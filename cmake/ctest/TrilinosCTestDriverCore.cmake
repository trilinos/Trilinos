
#
# Trilinos platform-independent ctest core driver
#
# This assumes that the script driving this will always be in the
# directory:
#
#   Trilinos/cmake/ctest
#
# which is set to CTEST_SCRIPT_DIRECTORY.
#
# All paths are relative to this directory.
#


CMAKE_MINIMUM_REQUIRED(VERSION 2.7.0 FATAL_ERROR)


# Get the base diretory for the Trilinos source.  We only assume that the
# CTest script that is being called is under Trilinos/cmake.
STRING(REGEX MATCH "(.+/Trilinos)/cmake" TRILINOS_CMAKE_DIR
  "${CTEST_SCRIPT_DIRECTORY}" )
IF("${TRILINOS_CMAKE_DIR}" STREQUAL "")
  STRING(REGEX MATCH "(.+)/cmake" TRILINOS_CMAKE_DIR
    "${CTEST_SCRIPT_DIRECTORY}" )
ENDIF()
  
MESSAGE("TRILINOS_CMAKE_DIR = ${TRILINOS_CMAKE_DIR}")

SET( CMAKE_MODULE_PATH
   "${TRILINOS_CMAKE_DIR}"
   "${TRILINOS_CMAKE_DIR}/utils"
   "${TRILINOS_CMAKE_DIR}/package_arch"
   )

#MESSAGE("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

INCLUDE(PrintVar)
INCLUDE(AssertDefined)
INCLUDE(PackageArchProcessPackagesAndDirsLists)


#
# Helper macros
#


FUNCTION(PRINT_VAR VAR)
  MESSAGE("${VAR} = '${${VAR}}'")
ENDFUNCTION()


MACRO(SET_DEFAULT VAR)
  IF ("${${VAR}}" STREQUAL "")
    SET(${VAR} ${ARGN})
  ENDIF()
ENDMACRO()


MACRO(SET_DEFAULT_AND_FROM_ENV VAR DEFAULT_VAL)

  SET_DEFAULT(${VAR} "${DEFAULT_VAL}")
  
  SET(ENV_${VAR} $ENV{${VAR}})
  IF (NOT "${ENV_${VAR}}" STREQUAL "")
    PRINT_VAR(ENV_${VAR})
    SET(${VAR} ${ENV_${VAR}})
  ENDIF()

  PRINT_VAR(${VAR})

ENDMACRO()



#
# Do some initial setup
#

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

# Get the host name
SITE_NAME(CTEST_SITE_DEFAULT)


#
# Select the list of packages
#
# OUTPUT: Sets Trilnos_DEFAULT_PACKAGES
#
# NOTE: This macro is used to cean up the main TRILINOS_CTEST_DRIVER()
# macro.
#

MACRO(SELECT_DEFAULT_TRILINOS_PACKAGES)

  # Must set this *before* reading the following file!
  SET(PROJECT_NAME Trilinos)

  INCLUDE(TrilinosPackages)
  #PRINT_VAR(Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)

  # Save the list of packages that might have been already sset
  SET( Trilinos_PACKAGES_SAVED ${Trilinos_PACKAGES})

  SET(${PROJECT_NAME}_ASSERT_MISSING_PACKAGES FALSE)
  SET(${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK TRUE)
  PACKAGE_ARCH_PROCESS_PACKAGES_AND_DIRS_LISTS()

  SET(Trilinos_PACKAGES_DEFAULT)

  FOREACH(PACKAGE ${Trilinos_PACKAGES})

    LIST(FIND Trilinos_EXCLUDE_PACKAGES ${PACKAGE} EXCLUDE_IDX)

    LIST(FIND Trilinos_ADDITIONAL_PACKAGES ${PACKAGE} ADDITIONAL_IDX)

    IF (EXCLUDE_IDX GREATER -1)
      # Exclude the package!
    ELSEIF (ADDITIONAL_IDX GREATER -1)
      LIST(APPEND Trilinos_PACKAGES_DEFAULT ${PACKAGE})
    ELSEIF (Trilinos_ENABLE_${PACKAGE} STREQUAL "")
      PACKAGE_ARCH_IMPLICIT_PACKAGE_ENABLE_IS_ALLOWED(
        ${PACKAGE} IMPLICIT_PACKAGE_ENABLE_ALLOWED)
      IF (IMPLICIT_PACKAGE_ENABLE_ALLOWED)
        LIST(APPEND Trilinos_PACKAGES_DEFAULT ${PACKAGE})
      ENDIF()
    ENDIF()

  ENDFOREACH()

  PRINT_VAR(Trilinos_PACKAGES_DEFAULT)

  # Reset the list of packages
  SET( Trilinos_PACKAGES ${Trilinos_PACKAGES_SAVED} )

ENDMACRO()


#
# This is the core extended ctest driver script code that is platform
# independent.  This script drives the testing process by doing an update and
# then configuring and building the packages one at a time.
#
# ToDo: Finish Documentation!
#

FUNCTION(TRILINOS_CTEST_DRIVER)

  #
  # Variables that can be set by the platform-specific code and reset
  # from the environment
  #
  
  # The type of test (e.g. Nightly, Experimental, Continuous)
  SET_DEFAULT_AND_FROM_ENV( CTEST_TEST_TYPE Nightly )
  
  # The default track to send the build to. This can be changed to send
  # the data to a different nightly grouping on the dashboard.
  SET_DEFAULT_AND_FROM_ENV(Trilinos_TRACK "")
  
  # The name of the site in the dashboard (almost never need to override this)
  SET_DEFAULT_AND_FROM_ENV( CTEST_SITE ${CTEST_SITE_DEFAULT} )

  # The root of the dasbhoard (almost never need to override this)
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

  # Select a default generator.  When the build tree is known and
  # exists, use its generator.
  SET(DEFAULT_GENERATOR "Unix Makefiles")
  IF(EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    FILE(STRINGS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" CACHE_CONTENTS)
    FOREACH(line ${CACHE_CONTENTS})
      IF("${line}" MATCHES "CMAKE_GENERATOR")
        STRING(REGEX REPLACE "(.*)=(.*)" "\\2" DEFAULT_GENERATOR "${line}")
      ENDIF()
    ENDFOREACH(line)
  ENDIF()
      
  SET_DEFAULT_AND_FROM_ENV( CTEST_CMAKE_GENERATOR ${DEFAULT_GENERATOR})

  # Do the CVS updates or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_UPDATES TRUE )
 
  # Generate the XML dependency output files or not
  SET_DEFAULT_AND_FROM_ENV( CTEST_GENERATE_DEPS_XML_OUTPUT_FILE FALSE )

  # Flags used on cvs when doing a CVS update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_ARGS "-q -z3")

  # Flags used on update when doing a CVS update
  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_OPTIONS "${Trilinos_BRANCH}")

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

  SET_DEFAULT_AND_FROM_ENV( Trilinos_ADDITIONAL_PACKAGES "" )

  SET_DEFAULT_AND_FROM_ENV( Trilinos_EXCLUDE_PACKAGES "" )

  SELECT_DEFAULT_TRILINOS_PACKAGES()

  SET_DEFAULT_AND_FROM_ENV( Trilinos_PACKAGES "${Trilinos_PACKAGES_DEFAULT}" )

  #
  # Setup and create the base dashboard directory if it is not created yet.
  #
  # NOTE: This is only used in general testing dashbaoard mode, not in local
  # testing mode.
  #

  SET( CTEST_SOURCE_NAME Trilinos )
  
  IF (CTEST_DASHBOARD_ROOT)
    SET( CTEST_BINARY_NAME BUILD )
    SET( CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
    SET( CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")
    IF (NOT EXISTS "${CTEST_DASHBOARD_ROOT}")
      MESSAGE("Creating the dashboard root directory \"${CTEST_DASHBOARD_ROOT}\" ...")
      FILE(MAKE_DIRECTORY "${CTEST_DASHBOARD_ROOT}")
    ENDIF()
  ENDIF()
  
  #
  # Some platform-independnet setup
  #
  
  INCLUDE("${TRILINOS_CMAKE_DIR}/../CTestConfig.cmake")
  SET(CMAKE_CACHE_CLEAN_FILE "${CTEST_BINARY_DIRECTORY}/CMakeCache.clean.txt")
  SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES};${CMAKE_CACHE_CLEAN_FILE}")
  SET(CTEST_USE_LAUNCHERS 1)
  
  #
  # Setup for the CVS update
  #

  IF (CTEST_DO_UPDATES)
    FIND_PACKAGE(CVS)
    #SET(CTEST_UPDATE_COMMAND "${CVS_EXECUTABLE} ${CTEST_UPDATE_ARGS}")
    SET(CTEST_UPDATE_COMMAND "${CVS_EXECUTABLE}")
    MESSAGE("CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")
    SET( CTEST_CHECKOUT_COMMAND
      "${CVS_EXECUTABLE} ${CTEST_UPDATE_ARGS} -d :ext:software.sandia.gov:/space/CVS co ${Trilinos_BRANCH} ${CTEST_SOURCE_NAME}" )
    MESSAGE("CTEST_CHECKOUT_COMMAND='${CTEST_CHECKOUT_COMMAND}'")
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
  
  #
  # Start up a new dashbaord
  #
  
  IF(Trilinos_TRACK)
    CTEST_START(${CTEST_TEST_TYPE} TRACK ${Trilinos_TRACK})
  ELSE()
    CTEST_START(${CTEST_TEST_TYPE})
  ENDIF()
  
  #
  # Do the VC update
  #

  IF (CTEST_DO_UPDATES)
    MESSAGE("Doing CVS update of '${CTEST_SOURCE_DIRECTORY}' ...")
    CTEST_UPDATE( SOURCE "${CTEST_SOURCE_DIRECTORY}"
      RETURN_VALUE  UPDATE_RETURN_VAL)
    MESSAGE("CTEST_UPDATE(...) returned '${UPDATE_RETURN_VAL}'")
  ENDIF()

  IF ("${UPDATE_RETURN_VAL}" LESS "0")
    MESSAGE("The VC update failed so submitting update and stopping ...") 
    IF (CTEST_DO_SUBMIT)
      CTEST_SUBMIT( PARTS update notes )
    ENDIF()
    RETURN()
  ENDIF()
  
  #
  # Tell CDash about the latest subproject dependencies:
  #
  
  PRINT_VAR(CTEST_DROP_SITE)
  IF (CTEST_DO_SUBMIT)
    CTEST_SUBMIT( FILES
      "${TRILINOS_CMAKE_DIR}/python/data/CDashSubprojectDependencies.xml"
      RETURN_VALUE SUBMIT_RETURN_VAL
      )
    MESSAGE("\nSubmitted subproject dependencies: Return='${SUBMIT_RETURN_VAL}'")
  ENDIF()
  
  #
  # loop over all Trilinos packages
  #
  
  MESSAGE("\nBegin incremental building and testing of Trilinos packages ...\n")
  
  SET(Trilinos_LAST_WORKING_PACKAGE)
  SET(Trilinos_FAILED_PACKAGES)
  
  FOREACH(PACKAGE ${Trilinos_PACKAGES})

    MESSAGE("\nKill all hanging Zoltan processes ...")
    EXECUTE_PROCESS(COMMAND killall -s 9 zdrive.exe)
  
    SET_PROPERTY(GLOBAL PROPERTY SubProject ${PACKAGE})
    SET_PROPERTY(GLOBAL PROPERTY Label ${PACKAGE})
  
    MESSAGE("\nCurrent Trilinos package: '${PACKAGE}'\n")
  
    #
    # Configure the package and its dependent packages
    #
  
    MESSAGE("\nConfiguring PACKAGE='${PACKAGE}'\n")
  
    # Create CONFIGURE_OPTIONS for this PACKAGE
    SET( CONFIGURE_OPTIONS
      "-DCTEST_USE_LAUNCHERS:BOOL=${CTEST_USE_LAUNCHERS}"
      "-DTrilinos_ENABLE_${PACKAGE}:BOOL=ON"
      "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
      "-DTrilinos_ENABLE_FRAMEWORK_UNIT_TESTS:BOOL=OFF"
      "-DTrilinos_ENABLE_TESTS:BOOL=ON"
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
    IF (CTEST_DO_COVERAGE_TESTING)
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_ENABLE_COVERAGE_TESTING:BOOL=ON")
    ENDIF()
    IF (BUILD_TYPE STREQUAL DEBUG)
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON")
    ENDIF()

    IF (DEFINED Trilinos_LAST_WORKING_PACKAGE)
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_ENABLE_${Trilinos_LAST_WORKING_PACKAGE}:BOOL=")
      SET(Trilinos_LAST_WORKING_PACKAGE)
    ENDIF()
    FOREACH(FAILED_PACKAGE ${Trilinos_FAILED_PACKAGES})
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
    ENDFOREACH()
    SET(CONFIGURE_OPTIONS ${CONFIGURE_OPTIONS}
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS} ${EXTRA_CONFIGURE_OPTIONS})
    MESSAGE("CONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")
  
    # Do the configure
    CTEST_CONFIGURE(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      OPTIONS "${CONFIGURE_OPTIONS}" # New option!
      RETURN_VALUE CONFIGURE_RETURN_VAL
      )
    MESSAGE("Generating the file CMakeCache.clean.txt ...")
    FILE(STRINGS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" CACHE_CONTENTS)
    message("CMAKE_CACHE_CLEAN_FILE = ${CMAKE_CACHE_CLEAN_FILE}")
    FILE(WRITE "${CMAKE_CACHE_CLEAN_FILE}")
    FOREACH(line ${CACHE_CONTENTS})
      # write lines that do not start with # or //
      IF(NOT "${line}" MATCHES "^(#|//)")
        FILE(APPEND "${CMAKE_CACHE_CLEAN_FILE}" "${line}\n")
      ENDIF()
    ENDFOREACH()
    # If the configure failed, add the package to the list
    # of failed packages
    IF (NOT "${CONFIGURE_RETURN_VAL}" EQUAL "0")
      MESSAGE("${PACKAGE} FAILED to configure")
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
    # If configure passed then try the build.  Otherwise, move on to
    # to the next package.
    #
  
    if ("${CONFIGURE_RETURN_VAL}" EQUAL "0")
  
      # Start by trying to build just the libraries for the current package
  
      SET( CTEST_BUILD_TARGET ${PACKAGE}_libs )
      MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")
      CTEST_BUILD (
        BUILD "${CTEST_BINARY_DIRECTORY}"
        RETURN_VALUE BUILD_LIBS_RETURN_VAL
        NUMBER_ERRORS BUILD_LIBS_NUM_ERRORS
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
  
        # Build the ALL target, but append the results to the last build.xml
        SET(CTEST_BUILD_TARGET)
        MESSAGE("\nBuild ALL target for '${PACKAGE}' ...\n")
        CTEST_BUILD(
          BUILD "${CTEST_BINARY_DIRECTORY}"
          NUMBER_ERRORS  BUILD_ALL_ERRORS
          RETURN_VALUE  BUILD_ALL_RETURN_VAL
          APPEND
          )
        MESSAGE("Build all: BUILD_ALL_ERRORS='${BUILD_ALL_ERRORS}',"
          "BUILD_ALL_RETURN_VAL='${BUILD_ALL_RETURN_VAL}'" )
  
        # Submit the build for all target
        IF (CTEST_DO_SUBMIT)
          CTEST_SUBMIT( PARTS build )
        ENDIF()
  
        IF (CTEST_DO_TEST)
          # Run the tests that match the ${PACKAGE} name 
          MESSAGE("\nRunning test for package '${PACKAGE}' ...\n")
          CTEST_TEST(BUILD "${CTEST_BINARY_DIRECTORY}"
            INCLUDE "^${PACKAGE}_"
            )
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
  
        # Remember this package so we can turn it off
        SET(Trilinos_LAST_WORKING_PACKAGE "${PACKAGE}")
  
      ELSE()
  
        MESSAGE("FAILED library build for package '${PACKAGE}'")
        LIST(APPEND Trilinos_FAILED_PACKAGES ${PACKAGE})
  
      ENDIF()
  
    ENDIF()

    IF (CTEST_DO_SUBMIT)
      MESSAGE("\nSubmit the update file that will trigger the notification email ...\n")
      CTEST_SUBMIT( PARTS update )
    ENDIF()

  ENDFOREACH(PACKAGE)
  
  IF(Trilinos_FAILED_PACKAGES)
    MESSAGE("\nFinal set of failed packages: '${Trilinos_FAILED_PACKAGES}'")
  ENDIF()

  MESSAGE("\nKill all hanging Zoltan processes ...")
  EXECUTE_PROCESS(COMMAND killall -s 9 zdrive.exe)
  
  MESSAGE("\nDone with the incremental building and testing of Trilinos packages!\n")

ENDFUNCTION()
