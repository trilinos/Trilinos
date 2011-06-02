
#
# Trilinos platform-independent ctest core driver
#
# This assumes that the outer CTest script driving this CTest code
# will always be in a sub-directory of:
#
#   Trilinos/cmake/ctest
#
# which is set to CTEST_SCRIPT_DIRECTORY automatically by CTest.
#
# All include and module paths are relative to this assumed directory
# structure.
#
# This script can be run from anywhere by by default and will find the
# right related CMake files to run but it requires that the client set
# CTEST_DASHBOARD_ROOT (or override this in the env) before running
# this script.  The varible CTEST_DASHBOARD_ROOT determines where the
# Trilinos source code will be cloned to and determines where the
# build directly will be put.
#

MESSAGE("")
MESSAGE("*******************************")
MESSAGE("*** TrilinosCTestDriverCore ***") 
MESSAGE("*******************************")
MESSAGE("")


CMAKE_MINIMUM_REQUIRED(VERSION 2.7.0 FATAL_ERROR)

# Must set this *before* reading the following file!
SET(PROJECT_NAME Trilinos)

# Get the base diretory for the Trilinos source.  We only assume that the
# CTest script that is being called is under Trilinos/cmake.
STRING(REGEX MATCH "(.+/Trilinos)/cmake" TRILINOS_CMAKE_DIR
  "${CTEST_SCRIPT_DIRECTORY}" )
IF("${TRILINOS_CMAKE_DIR}" STREQUAL "")
  STRING(REGEX MATCH "(.+)/cmake" TRILINOS_CMAKE_DIR
    "${CTEST_SCRIPT_DIRECTORY}" )
ENDIF()
  
MESSAGE("TRILINOS_CMAKE_DIR = ${TRILINOS_CMAKE_DIR}")

SET(TRILINOS_HOME_DIR "${TRILINOS_CMAKE_DIR}/..")
MESSAGE("TRILINOS_HOME_DIR = ${TRILINOS_HOME_DIR}")

SET( CMAKE_MODULE_PATH
   "${TRILINOS_CMAKE_DIR}"
   "${TRILINOS_CMAKE_DIR}/.."
   "${TRILINOS_CMAKE_DIR}/utils"
   "${TRILINOS_CMAKE_DIR}/package_arch"
   )

#MESSAGE("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

INCLUDE(PrintVar)
INCLUDE(SetDefaultAndFromEnv)
INCLUDE(AssertDefined)
INCLUDE(AppendSet)
INCLUDE(AppendStringVar)
INCLUDE(PackageArchGlobalMacros)

INCLUDE(TrilinosFindPythonInterp)
TRILINOS_FIND_PYTHON()
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

IF(WIN32)
  #Apparently FIND_PROGRAM looks for an exact match of the file name.
  #So even though "git clone ..." is valid to use on windows we need to give the
  #full name of the command we want to run.
  SET(GIT_NAME git.cmd)
ELSE()
  SET(GIT_NAME git)
ENDIF()
FIND_PROGRAM(GIT_EXE NAMES ${GIT_NAME})
MESSAGE("GIT_EXE=${GIT_EXE}")

IF(NOT GIT_EXE)
  QUEUE_ERROR("error: could not find git: GIT_EXE='${GIT_EXE}'")
ENDIF()
IF(NOT EXISTS "${GIT_EXE}")
  QUEUE_ERROR("error: GIT_EXE='${GIT_EXE}' does not exist")
ENDIF()

# Find svn

FIND_PROGRAM(SVN_EXE NAMES svn)
MESSAGE("SVN_EXE=${SVN_EXE}")
# If we don't find svn, no big deal unless we have an SVN repo!

# Get the host name

SITE_NAME(CTEST_SITE_DEFAULT)


########################
### Functions/Macros ###
########################


#
# Update or clone an extra repo
#

MACRO(CLONE_OR_UPDATE_EXTRAREPO  EXTRAREPO_NAME_IN  EXTRAREPO_DIR_IN
  EXTRAREPO_REPOTYPE_IN  EXTRAREPO_REPOURL_IN
  )

  #MESSAGE("CLONE_OR_UPDATE_EXTRAREPO: ${EXTRAREPO_NAME_IN} ${EXTRAREPO_REPOURL_IN}")

  SET(EXTRAREPO_SRC_DIR "${Trilinos_SOURCE_DIRECTORY}/${EXTRAREPO_DIR_IN}")
  SET(EXTRAREPO_UPDATE_OUT_FILE "${CTEST_BINARY_DIRECTORY}/${EXTRAREPO_NAME_IN}.update.out")
  #PRINT_VAR(EXTRAREPO_SRC_DIR)

  IF (NOT EXISTS "${EXTRAREPO_SRC_DIR}")

    MESSAGE("\n${EXTRAREPO_NAME_IN}: Doing initial ${EXTRAREPO_REPOTYPE_IN} clone/checkout from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_DIR_IN}' ...")
    IF (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
      SET(CMND_ARGS
        COMMAND "${GIT_EXE}" clone "${EXTRAREPO_REPOURL}" ${EXTRAREPO_DIR_IN}
        WORKING_DIRECTORY "${Trilinos_SOURCE_DIRECTORY}"
        OUTPUT_FILE "${EXTRAREPO_UPDATE_OUT_FILE}" )
    ELSEIF (${EXTRAREPO_REPOTYPE_IN} STREQUAL SVN)
      IF (NOT SVN_EXE)
        MESSAGE(SEND_ERROR "Error, could not find SVN executable!")
      ENDIF()
      SET(CMND_ARGS
        COMMAND "${SVN_EXE}" checkout "${EXTRAREPO_REPOURL}" ${EXTRAREPO_DIR_IN}
        WORKING_DIRECTORY "${Trilinos_SOURCE_DIRECTORY}"
        OUTPUT_FILE "${EXTRAREPO_UPDATE_OUT_FILE}" )
    ELSE()
      MESSAGE(SEND_ERROR "Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
    ENDIF()
    IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
      EXECUTE_PROCESS(${CMND_ARGS})
    ELSE()
      MESSAGE("EXECUTE_PROCESS(${CMND_ARGS})")
    ENDIF()

  ELSE()

    MESSAGE("\n${EXTRAREPO_NAME_IN}: Doing ${EXTRAREPO_REPOTYPE_IN} update from URL '${EXTRAREPO_REPOURL_IN}' to dir '${EXTRAREPO_SRC_DIR}' ...")
    IF (${EXTRAREPO_REPOTYPE_IN} STREQUAL GIT)
      SET(CMND_ARGS
        COMMAND "${GIT_EXE}" pull
        WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
        OUTPUT_FILE "${EXTRAREPO_UPDATE_OUT_FILE}" )
    ELSEIF (${EXTRAREPO_REPOTYPE_IN} STREQUAL SVN)
      SET(CMND_ARGS
        COMMAND "${SVN_EXE}" update
        WORKING_DIRECTORY "${EXTRAREPO_SRC_DIR}"
        OUTPUT_FILE "${EXTRAREPO_UPDATE_OUT_FILE}" )
    ELSE()
      MESSAGE(SEND_ERROR "Error, Invalid EXTRAREPO_REPOTYPE_IN='${EXTRAREPO_REPOTYPE_IN}'!")
    ENDIF()
    IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
      EXECUTE_PROCESS(${CMND_ARGS})
    ELSE()
      MESSAGE("EXECUTE_PROCESS(${CMND_ARGS})")
    ENDIF()

  ENDIF()

  # ToDo: Above, find a way to get the update info into the CDash submit to
  # be displayed on the cdash server.  If this is too hard to do then who
  # cares.  Eventually precopyright code will make it into the main
  # repository.

ENDMACRO()


#
# Select the set of extra Trilinos repositories
#

MACRO(SETUP_TRILINOS_EXTRAREPOS)

  MESSAGE("Reading the list of extra repos from ${${PROJECT_NAME}_EXTRAREPOS_FILE} ...")
  INCLUDE(${${PROJECT_NAME}_EXTRAREPOS_FILE})
  PACKAGE_ARCH_PROCESS_EXTRAREPOS_LISTS() # Sets ${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT
  SET_DEFAULT_AND_FROM_ENV(${PROJECT_NAME}_EXTRA_REPOSITORIES
    "${${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT}")
  #PRINT_VAR(${PROJECT_NAME}_EXTRA_REPOSITORIES)

ENDMACRO()


#
# Select the list of packages
#
# OUTPUT: Sets Trilinos_DEFAULT_PACKAGES
#
# NOTE: This macro is used to cean up the main TRILINOS_CTEST_DRIVER()
# macro.
#

MACRO(SETUP_TRILINOS_PACKAGES)

  # Here, we must point into the source tree of Trilinos just cloned (or
  # updated) and not the master Trilinos source dir tree for two reasons.
  # First, the list of core Trilinos packages may be more recent in what was
  # checked out.  Second, the extra repos do not even exist in the master
  # Trilinos source tree.
  IF (NOT ${PROJECT_NAME}_DEPS_HOME_DIR)
    SET(${PROJECT_NAME}_DEPS_HOME_DIR "${CTEST_SOURCE_DIRECTORY}")
  ENDIF()

  SET(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
  SET(${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK TRUE)
  SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES FALSE)
  SET(${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR "${CTEST_BINARY_DIRECTORY}")
  SET(PROJECT_HOME_DIR "${TRILINOS_HOME_DIR}")

  # Ignore missing extra repos in case someone messed up the git URL or
  # something.  We don't want a sloppy commit to bring down automated testing.
  SET(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES ON)

  PACKAGE_ARCH_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML()

  # When we get here, we will have the basic dependency structure set up
  # with only defaults set

  # Set this to "" so that it can be defined in ENABLE_MODIFIED_PACKAGES_ONLY()
  SET(Trilinos_ENABLE_ALL_PACKAGES "")

ENDMACRO()


#
# Select packages set by the input
#

MACRO(ENABLE_USER_SELECTED_PACKAGES)

  # 1) Set the enables for packages already set with Trilinos_PACKAGES_USER_SELECTED

  IF (NOT Trilinos_PACKAGES_USER_SELECTED)
    SET(Trilinos_ENABLE_ALL_PACKAGES ON)
  ELSE()
    FOREACH(PACKAGE ${Trilinos_PACKAGES_USER_SELECTED})
      MESSAGE("Enabling explicitly set package ${PACKAGE} ...")
      SET(Trilinos_ENABLE_${PACKAGE} ON)
    ENDFOREACH()
  ENDIF()

  # 2) Set extra package enables from Trilinos_ADDITIONAL_PACKAGES

  FOREACH(PACKAGE ${Trilinos_ADDITIONAL_PACKAGES})
    MESSAGE("Enabling explicitly set package ${PACKAGE} ...")
    SET(Trilinos_ENABLE_${PACKAGE} ON)
  ENDFOREACH()

ENDMACRO()


#
# Extract the list of changed files for the main Trilinos repo on put into
# an modified files file.
#

MACRO(GET_MODIFIED_TRILINOS_FILES  WORKING_DIR_IN  MODIFIED_FILES_FILE_NAME_IN)
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

MACRO(ENABLE_MODIFIED_PACKAGES_ONLY)

  #
  # A) Get the list of changed packages
  #

  SET(MODIFIED_FILES_FILE_NAME "${CTEST_BINARY_DIRECTORY}/modifiedFiles.txt")

  # A.1) Get changes from main Trilinos repo

  GET_MODIFIED_TRILINOS_FILES("${CTEST_SOURCE_DIRECTORY}" "${MODIFIED_FILES_FILE_NAME}")

  # A.2) Get changes from extra repos

  SET(EXTRAREPO_IDX 0)
  FOREACH(EXTRAREPO_NAME ${${PROJECT_NAME}_EXTRA_REPOSITORIES})

    LIST(GET Trilinos_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX} EXTRAREPO_DIR )
    LIST(GET Trilinos_EXTRA_REPOSITORIES_PACKSTATS ${EXTRAREPO_IDX} EXTRAREPO_PACKSTAT )
 
    # For now, only look for changes if it has packages.  Later, we need to
    # generalize this for the general extra repo case with deeper directory
    # and other than GIT (e.g. SVN with Dakota).  For example, we would like
    # to pick up changes to Dakota and therefore enable TriKota to build.
    IF (EXTRAREPO_PACKSTAT STREQUAL HASPACKAGES)

      SET(EXTRAREPO_SRC_DIR "${CTEST_SOURCE_DIRECTORY}/${EXTRAREPO_DIR}")
      SET(EXTRAREPO_MODIFIED_FILES_FILE_NAME
        "${CTEST_BINARY_DIRECTORY}/modifiedFiles.${EXTRAREPO_NAME}.txt")
  
      GET_MODIFIED_TRILINOS_FILES("${EXTRAREPO_SRC_DIR}" "${EXTRAREPO_MODIFIED_FILES_FILE_NAME}")
  
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
    MESSAGE(FATAL_ERROR "Error, Python must be enabled to map from modified files to packages!")
  ENDIF()

  EXECUTE_PROCESS(
    COMMAND ${PYTHON_EXECUTABLE}
      ${TRILINOS_CMAKE_DIR}/python/get-trilinos-packages-from-files-list.py
      --files-list-file=${MODIFIED_FILES_FILE_NAME}
      --deps-xml-file=${CTEST_BINARY_DIRECTORY}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}
    OUTPUT_VARIABLE MODIFIED_PACKAGES_LIST
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

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

  FOREACH(PACKAGE ${MODIFIED_PACKAGES_LIST})
    #PRINT_VAR(Trilinos_ENABLE_${PACKAGE})
    ASSERT_DEFINED(Trilinos_ENABLE_${PACKAGE})
    IF ("${Trilinos_ENABLE_${PACKAGE}}" STREQUAL "")
      MESSAGE("Enabling modified package: ${PACKAGE}")
      SET(Trilinos_ENABLE_${PACKAGE} ON)
    ELSE()
      MESSAGE("Not enabling explicitly disabled modified package: ${PACKAGE}")
    ENDIF()
  ENDFOREACH()

  FOREACH(PACKAGE ${FAILING_PACKAGES_LIST})
    IF ("${Trilinos_ENABLE_${PACKAGE}}" STREQUAL "")
      MESSAGE("Enabling previously failing package: ${PACKAGE}")
      SET(Trilinos_ENABLE_${PACKAGE} ON)
    ELSE()
      MESSAGE("Not enabling explicitly disabled previously failing package: ${PACKAGE}")
    ENDIF()
  ENDFOREACH()

  #
  # D) Print the final status
  #

  PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
    "\nDirectly modified or failing non-disabled packages that need to be tested" ON FALSE)

ENDMACRO()


#
# B) Exclude disabled packages from from Trilinos_EXCLUDE_PACKAGES
#
# NOTE: These disables need to dominate over the above enables so this code is
# after all the enable code has run
#

MACRO(DISABLE_EXCLUDED_PACKAGES)
  FOREACH(PACKAGE ${Trilinos_EXCLUDE_PACKAGES})
    MESSAGE("Disabling excluded package ${PACKAGE} ...")
    SET(Trilinos_ENABLE_${PACKAGE} OFF)
  ENDFOREACH()
ENDMACRO()


#
# Select the default generator.
#

MACRO(SELECT_DEFAULT_GENERATOR)
  # When the build tree is known and exists, use
  # its generator.
  SET(DEFAULT_GENERATOR "Unix Makefiles")
  IF(EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    FILE(STRINGS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" CACHE_CONTENTS)
    FOREACH(line ${CACHE_CONTENTS})
      IF("${line}" MATCHES "CMAKE_GENERATOR")
        STRING(REGEX REPLACE "(.*)=(.*)" "\\2" DEFAULT_GENERATOR "${line}")
      ENDIF()
    ENDFOREACH(line)
  ENDIF()
ENDMACRO()


#
# Call INITIALIZE_ERROR_QUEUE once at the top of TRILINOS_CTEST_DRIVER
#

MACRO(INITIALIZE_ERROR_QUEUE)
  SET(TRILINOS_CTEST_DRIVER_ERROR_QUEUE "")
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
  SET(TRILINOS_CTEST_DRIVER_ERROR_QUEUE
    ${TRILINOS_CTEST_DRIVER_ERROR_QUEUE} "${err_msg}")
ENDMACRO()


#
# Call REPORT_QUEUED_ERRORS once at the bottom of TRILINOS_CTEST_DRIVER
#

MACRO(REPORT_QUEUED_ERRORS)
  IF("${TRILINOS_CTEST_DRIVER_ERROR_QUEUE}" STREQUAL "")
    MESSAGE("TRILINOS_CTEST_DRIVER_ERROR_QUEUE is empty. All is well.")
  ELSE()
    MESSAGE("error: TRILINOS_CTEST_DRIVER_ERROR_QUEUE reports the following error message queue:")
    FOREACH(err_msg ${TRILINOS_CTEST_DRIVER_ERROR_QUEUE})
      MESSAGE("${err_msg}")
    ENDFOREACH()
  ENDIF()
ENDMACRO()


#
# Override CTEST_SUBMIT to detect failed submissions and track them as
# queued errors.
#

MACRO(CTEST_SUBMIT)
 
  # If using a recent enough ctest with RETRY_COUNT, use it to overcome
  # failed submits:
  #
  SET(retry_args "")
  IF("${CMAKE_VERSION}" VERSION_GREATER "2.8.2")
    SET(retry_args RETRY_COUNT 25 RETRY_DELAY 120)
    MESSAGE("info: using retry_args='${retry_args}' for _ctest_submit call")
  ENDIF()

  # Call the original CTEST_SUBMIT and pay attention to its RETURN_VALUE:
  #
  _CTEST_SUBMIT(${ARGN} ${retry_args} RETURN_VALUE rv)

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

#
# This is the core extended ctest driver script code that is platform
# independent.  This script drives the testing process by doing an update and
# then configuring and building the packages one at a time.
#
# ToDo: Finish Documentation!
#

FUNCTION(TRILINOS_CTEST_DRIVER)

  INITIALIZE_ERROR_QUEUE()

  SET( CTEST_SOURCE_NAME Trilinos )


  MESSAGE(
    "\n***\n"
    "\n*** Setting input options to default and reading from env ..."
    "\n***\n")
  
  # The type of test (e.g. Nightly, Experimental, Continuous)
  SET_DEFAULT_AND_FROM_ENV( CTEST_TEST_TYPE Nightly )
  
  # The default track to send the build to. This can be changed to send
  # the data to a different nightly grouping on the dashboard.
  # If the test type is set to Experimental though the track is forced
  # to "Experimental" this is so that we can have experimental tests 
  # on branches.
  SET_DEFAULT_AND_FROM_ENV(Trilinos_TRACK "")
  IF(CTEST_TEST_TYPE STREQUAL "Experimental" OR CTEST_TEST_TYPE STREQUAL "EXPERIMENTAL")
    SET(Trilinos_TRACK "Experimental")
  ENDIF()
 
  # The name of the site in the dashboard (almost never need to override this)
  SET_DEFAULT_AND_FROM_ENV( CTEST_SITE ${CTEST_SITE_DEFAULT} )

  # The root of the dasbhoard where Trilinos will be cloned and the
  # BUILD directory will be create (only override for separate testing)
  SET_DEFAULT_AND_FROM_ENV( CTEST_DASHBOARD_ROOT "" )

  # The build type (e.g. DEBUG, RELEASE, NONE)
  SET_DEFAULT_AND_FROM_ENV( BUILD_TYPE NONE )

  # Set the default compiler version
  SET_DEFAULT_AND_FROM_ENV(COMPILER_VERSION UNKNOWN)

  # The name of the build that appears in the dashbaord 
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
 
  # Generate the XML dependency output files or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_GENERATE_DEPS_XML_OUTPUT_FILE FALSE )

  # Flags used on git when doing a Git update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_ARGS "")

  # Flags used on update when doing a Git update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_OPTIONS "")

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
  # the test will be skipped) 
  SET_DEFAULT_AND_FROM_ENV( MPI_EXEC_MAX_NUMPROCS 4 )

  # How many tests ctest will spawn simultaneously
  SET_DEFAULT_AND_FROM_ENV( CTEST_PARALLEL_LEVEL 1 )

  # Turn off or change warnings-as-errors flag(s) (i.e. -Werror)
  SET_DEFAULT_AND_FROM_ENV( Trilinos_WARNINGS_AS_ERRORS_FLAGS "" )
  
  # Do coverage testing or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_COVERAGE_TESTING FALSE )

  # Command to run to get coverage results
  SET_DEFAULT_AND_FROM_ENV( CTEST_COVERAGE_COMMAND gcov )
  
  # Do memory testing (i.e. valgrind) or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_MEMORY_TESTING FALSE )

  # Command used to perform the memory testing (i.e. valgrind)
  SET_DEFAULT_AND_FROM_ENV( CTEST_MEMORYCHECK_COMMAND valgrind )
  
  # Submit the results to the dashboard or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_SUBMIT TRUE )

  SET_DEFAULT_AND_FROM_ENV( Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF )

  # List of additional packges that will be enabled over the current set
  # of all packagess (that would be set by Trilinos_ENABLE_ALL_PACKAGES).
  SET_DEFAULT_AND_FROM_ENV( Trilinos_ADDITIONAL_PACKAGES "" )

  # List of packages to not directly process.  NOTE: Listing these packages
  # here will *not* disable the package in the CMake build system.  To do
  # that, you will have to disable them in the variable
  # EXTRA_CONFIGURE_OPTIONS (set in your driver script.
  SET_DEFAULT_AND_FROM_ENV( Trilinos_EXCLUDE_PACKAGES "" )
  
  SET_DEFAULT_AND_FROM_ENV( Trilinos_BRANCH "" )

  IF(CTEST_TEST_TYPE STREQUAL "Nightly")
    SET_DEFAULT_AND_FROM_ENV( Trilinos_REPOSITORY_LOCATION "software.sandia.gov:/space/git/nightly/${CTEST_SOURCE_NAME}" )
  ELSE()
    SET_DEFAULT_AND_FROM_ENV( Trilinos_REPOSITORY_LOCATION "software.sandia.gov:/space/git/${CTEST_SOURCE_NAME}" )
  ENDIF()

  # Selct the Trilinos packages to enable (empty means to select all available)
  SET_DEFAULT_AND_FROM_ENV( Trilinos_PACKAGES "" )
  SET(Trilinos_PACKAGES_USER_SELECTED ${Trilinos_PACKAGES})
  SET(Trilinos_PACKAGES "")
  # Note: above, we have to keep the name Trilinos_PACKAGES to maintain
  # backward compatibility of this CTest script but we want to let
  # Trilinos_PACKAGES always be the full set of packages as defined by
  # the basic readin process

  # Override the location of the base directory where the package dependency
  # related files will be read relative to.  If left "", then this will be reset
  # to CTEST_SORUCE_DIRECTORY.
  SET_DEFAULT_AND_FROM_ENV(Trilinos_DEPS_HOME_DIR "")

  # Set the file that the extra repos will be read from
  #
  # NOTE: Here, we have no choice but to point into the master
  # Trilinos source treee because the local Trilinos sources have not
  # even been checked out yet!  Unless, of course, we are unit testing
  # in which case we will use whatever has been passed in.

  IF (NOT Trilinos_DEPS_HOME_DIR)
    SET(Trilinos_EXTRAREPOS_FILE_DEFAULT
      "${TRILINOS_CMAKE_DIR}/${Trilinos_EXTRA_EXTERNAL_REPOS_FILE_NAME}")
  ELSE()
    SET(Trilinos_EXTRAREPOS_FILE_DEFAULT
      "${Trilinos_DEPS_HOME_DIR}/cmake/${Trilinos_EXTRA_EXTERNAL_REPOS_FILE_NAME}")
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV(Trilinos_EXTRAREPOS_FILE ${Trilinos_EXTRAREPOS_FILE_DEFAULT})

  # Select the set of extra external repos to add in packages.
  # These are the same types as CTEST_TEST_TYPE (e.g. 'Continuous' and
  # 'Nightly').  This is set by default to ${CTEST_TEST_TYPE} can be
  # overridden independent of ${CTEST_TEST_TYPE} also.
  SET(Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT ${CTEST_TEST_TYPE})
  SET_DEFAULT_AND_FROM_ENV( Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE
     "${Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_DEFAULT}" )

  # Set as part of CI testing in order to only enable modified packages
  SET_DEFAULT_AND_FROM_ENV( CTEST_ENABLE_MODIFIED_PACKAGES_ONLY OFF )

  # Set if implicitly enabled packages should be explicitly handled
  IF (CTEST_ENABLE_MODIFIED_PACKAGES_ONLY AND NOT CTEST_START_WITH_EMPTY_BINARY_DIRECTORY)
    SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT FALSE )
  ELSE()
    SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT TRUE )
  ENDIF()
  SET_DEFAULT_AND_FROM_ENV( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES
    ${CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES_DEFAULT})

  MESSAGE(
    "\n***"
    "\n*** Setting unit testing input options to default and reading from env ..."
    "\n***\n")

  SET_DEFAULT_AND_FROM_ENV( CTEST_DEPENDENCY_HANDLING_UNIT_TESTING FALSE )

  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_RETURN_VAL 0 )

  IF (CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    SET(GIT_EXE /somebasedir/git)
    SET(SVN_EXE /someotherbasedir/svn)
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

  # Set override hook for unit testing
  SET_DEFAULT_AND_FROM_ENV( Trilinos_SOURCE_DIRECTORY ${CTEST_SOURCE_DIRECTORY} )

  # Must be set here after CTEST_BINARY_DIRECTORY is set!
  SET(FAILED_PACKAGES_FILE_NAME "${CTEST_BINARY_DIRECTORY}/failedPackages.txt")

  
  #
  # Some platform-independent setup
  #
  
  INCLUDE("${TRILINOS_CMAKE_DIR}/../CTestConfig.cmake")
  SET(CMAKE_CACHE_CLEAN_FILE "${CTEST_BINARY_DIRECTORY}/CMakeCache.clean.txt")
  SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES};${CMAKE_CACHE_CLEAN_FILE}")
  SET(CTEST_USE_LAUNCHERS 1)
  
  #
  # Setup for the VC update
  #

  IF (CTEST_DO_UPDATES)

    SET(UPDATE_TYPE "git")
    MESSAGE("UPDATE_TYPE = '${UPDATE_TYPE}'")

    SET(CTEST_UPDATE_COMMAND "${GIT_EXE}")
    MESSAGE("CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")
    IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
      MESSAGE("${CTEST_SOURCE_DIRECTORY} does not exist so setting up for an initial checkout")
      SET( CTEST_CHECKOUT_COMMAND
        "\"${GIT_EXE}\" clone ${CTEST_UPDATE_ARGS} ${Trilinos_REPOSITORY_LOCATION}" )
      MESSAGE("CTEST_CHECKOUT_COMMAND='${CTEST_CHECKOUT_COMMAND}'")
    ELSE()
      MESSAGE("${CTEST_SOURCE_DIRECTORY} exists so skipping the initial checkout.")
    ENDIF()
  ENDIF() 


  #
  # Empty out the binary directory
  #

  IF (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY)
    MESSAGE("Cleaning out binary directory '${CTEST_BINARY_DIRECTORY}' ...")
    CTEST_EMPTY_BINARY_DIRECTORY("${CTEST_BINARY_DIRECTORY}")
  ELSEIF (CTEST_WIPE_CACHE)
    SET(CACHE_FILE_NAME "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    IF (EXISTS "${CACHE_FILE_NAME}")
      MESSAGE("Removing existing cache file '${CACHE_FILE_NAME}' ...")
      FILE(REMOVE "${CACHE_FILE_NAME}")
    ENDIF()
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Read in the set of extra repos ..."
    "\n***\n")

  SETUP_TRILINOS_EXTRAREPOS()

  # NOTE: You have to set up the set of extra repos before you can read the
  # Dependencies.cmake files since the extra repos must be cloned first.


  MESSAGE(
    "\n***"
    "\n*** Start up a new dashboard ..."
    "\n***\n")

  PRINT_VAR(CTEST_TEST_TYPE)
  PRINT_VAR(Trilinos_TRACK)

  IF(Trilinos_TRACK)
    CTEST_START(${CTEST_TEST_TYPE} TRACK ${Trilinos_TRACK})
  ELSE()
    CTEST_START(${CTEST_TEST_TYPE})
  ENDIF()

  
  MESSAGE(
    "\n***"
    "\n*** Update the source code repositories ..."
    "\n***\n")

  SET(UPDATE_FAILED TRUE)

  IF (CTEST_DO_UPDATES)

    IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
      QUEUE_ERROR("error: source directory does not exist just prior to CTEST_UPDATE call -- initial checkout did not work")
      REPORT_QUEUED_ERRORS()
      RETURN()
    ENDIF()

    MESSAGE("\nDoing GIT update of '${CTEST_SOURCE_DIRECTORY}' ...")
    CTEST_UPDATE_WRAPPER( SOURCE "${CTEST_SOURCE_DIRECTORY}"
      RETURN_VALUE  UPDATE_RETURN_VAL)
    MESSAGE("CTEST_UPDATE(...) returned '${UPDATE_RETURN_VAL}'")

    SET(EXTRAREPO_IDX 0)
    FOREACH(EXTRAREPO_NAME ${Trilinos_EXTRA_REPOSITORIES})
      LIST(GET Trilinos_EXTRA_REPOSITORIES_DIRS ${EXTRAREPO_IDX} EXTRAREPO_DIR )
      LIST(GET Trilinos_EXTRA_REPOSITORIES_REPOTYPES ${EXTRAREPO_IDX} EXTRAREPO_REPOTYPE )
      LIST(GET Trilinos_EXTRA_REPOSITORIES_REPOURLS ${EXTRAREPO_IDX} EXTRAREPO_REPOURL )
      CLONE_OR_UPDATE_EXTRAREPO(${EXTRAREPO_NAME} ${EXTRAREPO_DIR}
        ${EXTRAREPO_REPOTYPE} ${EXTRAREPO_REPOURL})
      MATH(EXPR EXTRAREPO_IDX "${EXTRAREPO_IDX}+1")
    ENDFOREACH()

    # Setting branch switch to success in case we are not doing a switch to a
    # different branch.

    SET(GIT_CHECKOUT_RETURN_VAL "0")

    IF(Trilinos_BRANCH AND NOT "${UPDATE_RETURN_VAL}" LESS "0")

      MESSAGE("Doing switch to branch ${Trilinos_BRANCH}")

      EXECUTE_PROCESS(COMMAND ${GIT_EXE} checkout ${Trilinos_BRANCH}
        WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
        RESULT_VARIABLE GIT_CHECKOUT_RETURN_VAL
        OUTPUT_VARIABLE BRANCH_OUTPUT
        ERROR_VARIABLE  BRANCH_ERROR
      )

      IF(NOT "${GIT_CHECKOUT_RETURN_VAL}" EQUAL "0")
        EXECUTE_PROCESS(COMMAND ${GIT_EXE} checkout --track origin/${Trilinos_BRANCH}
          WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
          RESULT_VARIABLE GIT_CHECKOUT_RETURN_VAL
          OUTPUT_VARIABLE BRANCH_OUTPUT
          ERROR_VARIABLE  BRANCH_ERROR
        )
      ENDIF()

      IF(NOT "${GIT_CHECKOUT_RETURN_VAL}" EQUAL "0")
        MESSAGE("Switch to branch ${Trilinos_BRANCH} failed with error code ${GIT_CHECKOUT_RETURN_VAL}")
        QUEUE_ERROR("Switch to branch ${Trilinos_BRANCH} failed with error code ${GIT_CHECKOUT_RETURN_VAL}")
      ENDIF()
      #Apparently the successful branch switch is also written to stderr.
      MESSAGE("${BRANCH_ERROR}")

    ENDIF()

    IF ("${UPDATE_RETURN_VAL}" LESS "0" OR NOT "${GIT_CHECKOUT_RETURN_VAL}" EQUAL "0")
      SET(UPDATE_FAILED TRUE)
    ELSE()
      SET(UPDATE_FAILED FALSE)
    ENDIF()

  ELSE()

     MESSAGE("Skipping the update by request!")

     SET(UPDATE_FAILED FALSE)

  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Read in the set of packages and their dependencies ..."
    "\n***\n")

  # NOTE: You must read the Dependencies.cmake files *after* you have cloned
  # (or updated) all of the code!

  SETUP_TRILINOS_PACKAGES()

  SET(CDASH_SUBPROJECT_XML_FILE
    "${CTEST_BINARY_DIRECTORY}/${Trilinos_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}")
  PRINT_VAR(CDASH_SUBPROJECT_XML_FILE)

  DISABLE_EXCLUDED_PACKAGES()

  IF (NOT CTEST_ENABLE_MODIFIED_PACKAGES_ONLY)
    MESSAGE(
      "\n***"
      "\n*** Determining what packages to enable based what was set in the input ..."
      "\n***\n")
    ENABLE_USER_SELECTED_PACKAGES()
  ELSE()
    MESSAGE(
      "\n***"
      "\n*** Determining what packages to enable based on what changed ..."
      "\n***\n")
    ENABLE_MODIFIED_PACKAGES_ONLY()
    SET(Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES ON)
  ENDIF()


  MESSAGE(
    "\n***"
    "\n*** Adjust the package dependencies to enable upstream and (optionally) downstream packages ..."
    "\n***"
    )

  SET(Trilinos_ENABLE_TESTS ON)
  SET(Trilinos_ENABLE_EXAMPLES ON)
  SET(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES ON)
  SET(DO_PROCESS_MPI_ENABLES FALSE) # Should not be needed but CMake is messing up
  PACKAGE_ARCH_ADJUST_AND_PRINT_PACKAGE_DEPENDENCIES() # Sets Trilinos_NUM_ENABLED_PACKAGES

  
  MESSAGE(
    "\n***"
    "\n*** Determine if to go ahead with configure, build, test ..."
    "\n***")

  IF (CTEST_ENABLE_MODIFIED_PACKAGES_ONLY)
    IF (MODIFIED_PACKAGES_LIST)
      MESSAGE("\nMODIFIED_PACKAGES_LIST='${MODIFIED_PACKAGES_LIST}'"
        ":  Found modified packages, processing enabled packages!\n")
    ELSE()
      MESSAGE("\nMODIFIED_PACKAGES_LIST='${MODIFIED_PACKAGES_LIST}'"
        ":  No modified packages to justify continuous integration test iteration!\n")
      REPORT_QUEUED_ERRORS()
      RETURN()
    ENDIF()
  ELSE()
    MESSAGE(
      "\nCTEST_ENABLE_MODIFIED_PACKAGES_ONLY=${CTEST_ENABLE_MODIFIED_PACKAGES_ONLY}"
      "  Running in regular mode, processing all enabled packages!\n")
  ENDIF()

  IF (Trilinos_NUM_ENABLED_PACKAGES GREATER 0)
    MESSAGE(
      "\nTrilinos_NUM_ENABLED_PACKAGES=${Trilinos_NUM_ENABLED_PACKAGES}:"
      "  Configuring packages!\n")
  ELSE()
    MESSAGE(
      "\nTrilinos_NUM_ENABLED_PACKAGES=${Trilinos_NUM_ENABLED_PACKAGES}:"
      "  Exiting the script!\n")
    REPORT_QUEUED_ERRORS()
    RETURN()
  ENDIF()

  
  MESSAGE(
    "\n***"
    "\n*** Uploading update, notes, and the subproject dependencies XML files ..."
    "\n***\n"
    )

  # Note: We must only do the submit after we have decided if there are any
  # packages to enable or not and otherwise exit the script!

  IF (UPDATE_FAILED)
    MESSAGE("The VC update failed so submitting update and stopping ...") 
    IF (CTEST_DO_SUBMIT)
      CTEST_SUBMIT( PARTS update notes )
    ENDIF()
    REPORT_QUEUED_ERRORS()
    RETURN()
  ENDIF()

  IF (CTEST_DO_SUBMIT)
    CTEST_SUBMIT( FILES ${CDASH_SUBPROJECT_XML_FILE})
    MESSAGE("\nSubmitted subproject dependencies XML file!")
  ELSE()
    MESSAGE("\nSkipping submitted subproject dependencies XML file on request!")
  ENDIF()

  
  MESSAGE(
    "\n***"
    "\n*** Loop through Trilinos packages to configure, build, and test ..."
    "\n***")
  
  SET(Trilinos_LAST_CONFIGURED_PACKAGE)
  SET(Trilinos_FAILED_LIB_BUILD_PACKAGES)
  SET(Trilinos_FAILED_PACKAGES)
  SET(PACKAGE_IDX 0)
  
  FOREACH(PACKAGE ${Trilinos_PACKAGES})

    MESSAGE("")
    MESSAGE("${PACKAGE_IDX}) Procesing current package ${PACKAGE}: libs='${Trilinos_ENABLE_${PACKAGE}}', tests='${${PACKAGE}_ENABLE_TESTS}'")
    MESSAGE("")

    IF (Trilinos_ENABLE_${PACKAGE} AND NOT ${PACKAGE}_ENABLE_TESTS AND
     NOT CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES
     )

      MESSAGE("Not enabling implicitly enabled package ${PACKAGE} on request!")

    ELSEIF (Trilinos_ENABLE_${PACKAGE})

      SET_PROPERTY(GLOBAL PROPERTY SubProject ${PACKAGE})
      SET_PROPERTY(GLOBAL PROPERTY Label ${PACKAGE})
   
      #
      # A) Configure the package and its dependent packages
      #
    
      MESSAGE("Configuring PACKAGE='${PACKAGE}'")
    
      # Create CONFIGURE_OPTIONS for this PACKAGE
      SET( CONFIGURE_OPTIONS
        "-DCTEST_USE_LAUNCHERS:BOOL=${CTEST_USE_LAUNCHERS}"
        "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
        "-DTrilinos_ENABLE_TESTS:BOOL=${${PACKAGE}_ENABLE_TESTS}"
        "-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=${Trilinos_WARNINGS_AS_ERRORS_FLAGS}"
        "-DTrilinos_ALLOW_NO_PACKAGES:BOOL=OFF"
        )
      IF (NOT CTEST_GENERATE_DEPS_XML_OUTPUT_FILE)
        LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_DEPS_XML_OUTPUT_FILE:FILEPATH=")
      ENDIF()
      IF (Trilinos_ENABLE_SECONDARY_STABLE_CODE)
        LIST(APPEND CONFIGURE_OPTIONS
          "-DTrilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON")
      ENDIF()
      IF (MPI_EXEC_MAX_NUMPROCS)
        LIST(APPEND CONFIGURE_OPTIONS
          "-DMPI_EXEC_MAX_NUMPROCS:STRING=${MPI_EXEC_MAX_NUMPROCS}")
      ENDIF()
      IF (CTEST_DO_COVERAGE_TESTING)
        LIST(APPEND CONFIGURE_OPTIONS
          "-DTrilinos_ENABLE_COVERAGE_TESTING:BOOL=ON")
      ENDIF()
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_EXTRAREPOS_FILE:STRING=${Trilinos_EXTRAREPOS_FILE}")
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=${Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE}")
      IF (DEFINED Trilinos_LAST_CONFIGURED_PACKAGE)
        LIST(APPEND CONFIGURE_OPTIONS
          "-DTrilinos_ENABLE_${Trilinos_LAST_CONFIGURED_PACKAGE}:BOOL=")
        SET(Trilinos_LAST_CONFIGURED_PACKAGE)
      ENDIF()
      FOREACH(FAILED_PACKAGE ${Trilinos_FAILED_LIB_BUILD_PACKAGES})
        LIST(APPEND CONFIGURE_OPTIONS
          "-DTrilinos_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
      ENDFOREACH()
      SET(CONFIGURE_OPTIONS ${CONFIGURE_OPTIONS}
        ${EXTRA_SYSTEM_CONFIGURE_OPTIONS} ${EXTRA_CONFIGURE_OPTIONS})
      LIST(APPEND CONFIGURE_OPTIONS # Package enable must be at the very end to override other stuff!
         "-DTrilinos_ENABLE_${PACKAGE}:BOOL=ON" )
      MESSAGE("\nCONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")

      # Remember this package so we can set its enable to "" next time
      SET(Trilinos_LAST_CONFIGURED_PACKAGE "${PACKAGE}")

      #
      # B) Configure the package and its dependent packages
      #

      IF (NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    
        CTEST_CONFIGURE(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          OPTIONS "${CONFIGURE_OPTIONS}" # New option!
          RETURN_VALUE CONFIGURE_RETURN_VAL
          )
    
        MESSAGE("Generating the file CMakeCache.clean.txt ...")
        FILE(STRINGS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" CACHE_CONTENTS)
        MESSAGE("CMAKE_CACHE_CLEAN_FILE = ${CMAKE_CACHE_CLEAN_FILE}")
        SET(CMAKE_CACHE_CLEAN_FILE_STR "")
        FOREACH(line ${CACHE_CONTENTS})
          # write lines that do not start with # or //
          IF(NOT "${line}" MATCHES "^(#|//)")
            APPEND_STRING_VAR(CMAKE_CACHE_CLEAN_FILE_STR "${line}\n")
          ENDIF()
        ENDFOREACH()
        FILE(WRITE "${CMAKE_CACHE_CLEAN_FILE}" ${CMAKE_CACHE_CLEAN_FILE_STR})
    
        # If the configure failed, add the package to the list
        # of failed packages
        IF (NOT "${CONFIGURE_RETURN_VAL}" EQUAL "0")
          MESSAGE("${PACKAGE} FAILED to configure")
          LIST(APPEND Trilinos_FAILED_LIB_BUILD_PACKAGES ${PACKAGE})
          LIST(APPEND Trilinos_FAILED_PACKAGES ${PACKAGE})
        ELSE()
          # load target properties and test keywords
          CTEST_READ_CUSTOM_FILES(BUILD "${CTEST_BINARY_DIRECTORY}")
          # Overridde from this file!
          INCLUDE("${TRILINOS_CMAKE_DIR}/../CTestConfig.cmake")
        ENDIF()
      
        # Submit configure results and the notes to the dashboard 
        IF (CTEST_DO_SUBMIT)
          MESSAGE("\nSubmitting configure and notes ...")
          CTEST_SUBMIT( PARTS configure notes )
        ENDIF()
      
        #
        # C) If configure passed then try the build.  Otherwise, move on to
        # to the next package.
        #

      ENDIF()
    
      IF ("${CONFIGURE_RETURN_VAL}" EQUAL "0" AND NOT CTEST_DEPENDENCY_HANDLING_UNIT_TESTING)
    
        # Start by trying to build just the libraries for the current package
    
        SET( CTEST_BUILD_TARGET ${PACKAGE}_libs )
        MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")
        CTEST_BUILD(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          RETURN_VALUE  BUILD_LIBS_RETURN_VAL
          NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
          APPEND
          )
        MESSAGE("Build return: RETURN_VALUE=${BUILD_LIBS_RETURN_VAL},"
          " NUMBER_ERRORS=${BUILD_LIBS_NUM_ERRORS}")
    
        # Determine if the build failed or not.
    
        SET(BUILD_LIBS_SUCCESS FALSE)
        IF ("${BUILD_LIBS_NUM_ERRORS}" EQUAL "0" AND "${BUILD_LIBS_RETURN_VAL}" EQUAL "0")
          SET(BUILD_LIBS_SUCCESS TRUE)
        ENDIF()
        # Above: Since make -i is used BUILD_LIBS_RETURN_VAL might be 0, but
        # if there are errors the build should fail, so
        # both BUILD_LIBS_RETURN_VAL and BUILD_LIBS_NUM_ERRORS should be 0 for a good build
        # and for the all target to be built.
    
        # Submit the library build results to the dashboard
    
        IF (CTEST_DO_SUBMIT)
          CTEST_SUBMIT( PARTS build )
        ENDIF()
    
        # If the build of the libraries passed, then go on the build
        # the tests/examples and run them.

        IF (BUILD_LIBS_SUCCESS)
  
          SET(BUILD_OR_TEST_FAILED FALSE)
    
          # Build the ALL target, but append the results to the last build.xml
          SET(CTEST_BUILD_TARGET)
          MESSAGE("\nBuild ALL target for '${PACKAGE}' ...\n")
          CTEST_BUILD(
            BUILD "${CTEST_BINARY_DIRECTORY}"
            RETURN_VALUE  BUILD_ALL_RETURN_VAL
            NUMBER_ERRORS  BUILD_ALL_NUM_ERRORS
            APPEND
            )
          MESSAGE("Build all: BUILD_ALL_NUM_ERRORS='${BUILD_ALL_NUM_ERRORS}',"
            "BUILD_ALL_RETURN_VAL='${BUILD_ALL_RETURN_VAL}'" )
    
          IF (NOT "${BUILD_LIBS_NUM_ERRORS}" EQUAL "0" OR NOT "${BUILD_LIBS_RETURN_VAL}" EQUAL "0")
            SET(BUILD_OR_TEST_FAILED TRUE)
          ENDIF()
  
          # Submit the build for all target
          IF (CTEST_DO_SUBMIT)
            CTEST_SUBMIT( PARTS build )
          ENDIF()
    
          IF (CTEST_DO_TEST)
            # Remove the LastTestsFailed log so we can detect if there are any
            # failed tests.
            SET(TEST_TMP_DIR "${CTEST_BINARY_DIRECTORY}/Testing/Temporary")
            FILE(GLOB logfiles "${TEST_TMP_DIR}/LastTestsFailed*.log")
            FOREACH(logfile ${logfiles})
              FILE(REMOVE "${logfile}")
            ENDFOREACH()
            # Run the tests that match the ${PACKAGE} name 
            MESSAGE("\nRunning test for package '${PACKAGE}' ...\n")
            CTEST_TEST(
              BUILD "${CTEST_BINARY_DIRECTORY}"
              PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
              INCLUDE "^${PACKAGE}_"
              #NUMBER_FAILED  TEST_NUM_FAILED
              )
            # See if a 'LastTestsFailed*.log' file exists to determine if there
            # are failed tests
            FILE(GLOB FAILED_TEST_LOG_FILE "${TEST_TMP_DIR}/LastTestsFailed*.log")
            PRINT_VAR(FAILED_TEST_LOG_FILE)
            IF (FAILED_TEST_LOG_FILE)
              SET(BUILD_OR_TEST_FAILED TRUE)
            ENDIF()
            # 2009/12/05: ToDo: We need to add an argument to CTEST_TEST(...) 
            # called something like 'NUMBER_FAILED numFailedTests' to allow us
            # to detect when the tests have filed.
            #IF (TEST_NUM_FAILED GREATER 0)
            #  SET(BUILD_OR_TEST_FAILED TRUE)
            #ENDIF()
            IF (CTEST_DO_SUBMIT)
              CTEST_SUBMIT( PARTS Test )
            ENDIF()
          ENDIF()
  
          IF (CTEST_DO_COVERAGE_TESTING)
            MESSAGE("\nRunning coverage for package '${PACKAGE}' ...\n")
            CTEST_COVERAGE(
              BUILD "${CTEST_BINARY_DIRECTORY}"
              LABELS ${PACKAGE}
              )
            IF (CTEST_DO_SUBMIT)
              CTEST_SUBMIT( PARTS Coverage )
            ENDIF()
          ENDIF() 
   
          IF (CTEST_DO_MEMORY_TESTING)
            MESSAGE("\nRunning memory testing for package '${PACKAGE}' ...\n")
            CTEST_MEMCHECK(BUILD "${CTEST_BINARY_DIRECTORY}")
            IF (CTEST_DO_SUBMIT)
              CTEST_SUBMIT( PARTS Memcheck )
            ENDIF()
          ENDIF()
  
          IF (BUILD_OR_TEST_FAILED)
            LIST(APPEND Trilinos_FAILED_PACKAGES ${PACKAGE})
          ENDIF()
    
        ELSE()
    
          MESSAGE("FAILED library build for package '${PACKAGE}'")
          LIST(APPEND Trilinos_FAILED_LIB_BUILD_PACKAGES ${PACKAGE})
          LIST(APPEND Trilinos_FAILED_PACKAGES ${PACKAGE})
    
        ENDIF()
    
      ENDIF()
  
      IF (CTEST_DO_SUBMIT)
        MESSAGE("\nSubmit the update file that will trigger the notification email ...\n")
        CTEST_SUBMIT( PARTS update )
      ENDIF()

    ELSE()

      MESSAGE("Package ${PACKAGE} is disabled, skipping configure, build, test ...")

    ENDIF()

    MATH(EXPR PACKAGE_IDX "${PACKAGE_IDX}+1")

  ENDFOREACH(PACKAGE)
  
  IF(Trilinos_FAILED_LIB_BUILD_PACKAGES)
    MESSAGE(
      "\nFinal set packages that failed to configure or have the libraries build:"
      " '${Trilinos_FAILED_LIB_BUILD_PACKAGES}'")
  ENDIF()

  IF(Trilinos_FAILED_PACKAGES)
    MESSAGE(
      "\nFinal set packages that had any failures: '${Trilinos_FAILED_PACKAGES}'")
  ENDIF()

  # Write a file listing the packages that failed.  This will be read in on the next CI
  # iteration since these packages must be enabled
  FILE(WRITE "${FAILED_PACKAGES_FILE_NAME}" "${Trilinos_FAILED_PACKAGES}\n")

  # This is no longer necessary with CMake 2.8.1
  #MESSAGE("\nKill all hanging Zoltan processes ...")
  #EXECUTE_PROCESS(COMMAND killall -s 9 zdrive.exe)
  
  MESSAGE("\nDone with the incremental building and testing of Trilinos packages!\n")

  REPORT_QUEUED_ERRORS()

ENDFUNCTION()
