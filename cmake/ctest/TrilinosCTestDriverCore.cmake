
#
# Helper macros
#


FUNCTION(ASSERT_DEFINED VARIBLE_NAME)
  IF(NOT DEFINED ${VARIBLE_NAME})
    MESSAGE(SEND_ERROR "Error, the variable ${VARIBLE_NAME} is not defined!")
  ENDIF()
ENDFUNCTION()


FUNCTION(PRINT_VAR VAR)
  MESSAGE("${VAR} = '${${VAR}}'")
ENDFUNCTION()


MACRO(SET_DEFAULT VAR DEFAULT_VAL)
  
  IF ("${${VAR}}" STREQUAL "")
    SET(${VAR} ${DEFAULT_VAL})
  ENDIF()

ENDMACRO()


MACRO(SET_DEFAULT_AND_FROM_ENV VAR DEFAULT_VAL)

  SET_DEFAULT(${VAR} ${DEFAULT_VAL})
  
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
FIND_PROGRAM(UNAME_EXE NAMES uname)
EXECUTE_PROCESS(
  COMMAND ${UNAME_EXE}
  OUTPUT_VARIABLE HOST_TYPE
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# Get the host name

SITE_NAME(CTEST_SITE_DEFAULT)


#
# This is the core extended ctest driver script code that is platform
# independent.  This script drives the testing process by doing an update and
# then configuring and building the packages one at a time.
#
# ToDo: Finish Documentation!
#

MACRO(TRILINOS_CTEST_DRIVER)

  #
  # Variables that can be set by the platform-specific code and reset
  # from the environment
  #
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_TEST_TYPE Nightly )
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_SITE ${CTEST_SITE_DEFAULT} )
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE )
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_UPDATES TRUE )

  SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_FLAGS "-j2")

  SET_DEFAULT_AND_FROM_ENV( CTEST_UPDATE_ARGS "-q -z3")
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_BUILD TRUE )
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_TEST TRUE )
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_COVERAGE_TESTING FALSE )
  
  SET_DEFAULT_AND_FROM_ENV( CTEST_DO_MEMORY_TESTING FALSE )
  
  #
  # Some Platform-independnet setup
  #
  
  INCLUDE(${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)
  SET(CTEST_USE_LAUNCHERS 1)
  SET(CTEST_NOTES_FILES "${CTEST_NOTES_FILES};${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
  
  #
  # Setup for the CVS update
  #

  IF (CTEST_DO_UPDATES)
    FIND_PACKAGE(CVS)
    SET(CTEST_UPDATE_COMMAND "${CVS_EXECUTABLE} ${CTEST_UPDATE_ARGS}")
    MESSAGE("CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")
  ENDIF() 
 
  # 2009/02/18: rabartl: Above: If this update fails, how will I figure this
  # out?  Will the dashboard get updated?
  
  #
  # Get the list of Trilinos packages to do the tests on
  #
  
  SET( Trilinos_PACKAGES_DEFAULT
    Teuchos
    RTOp
    Epetra
    Thyra
    )
  # ToDo: Read this list from TrilinosPackages.cmake

  SET_DEFAULT_AND_FROM_ENV( Trilinos_PACKAGES "${Trilinos_PACKAGES_DEFAULT}" )
  
  #
  # Empty out the binary directory
  #
  
  IF (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY)
    MESSAGE("Cleaning out binary directory '${CTEST_BINARY_DIRECTORY}'")
    CTEST_EMPTY_BINARY_DIRECTORY("${CTEST_BINARY_DIRECTORY}")
  ENDIF()
  
  #
  # Start up a new dashbaord
  #
  
  CTEST_START(${CTEST_TEST_TYPE})

  #
  # Do the VC update
  #

  IF (CTEST_DO_UPDATES)
    MESSAGE("Doing CVS update of '${CTEST_SOURCE_DIRECTORY}' ...")
    CTEST_UPDATE( SOURCE "${CTEST_SOURCE_DIRECTORY}"
      RETURN_VALUE  UPDATE_RETURN_VAL)
    MESSAGE("Updated return = '${UPDATE_RETURN_VAL}'")
  ENDIF()
  
  #
  # Tell CDash about the latest subproject dependencies:
  #
  
  CTEST_SUBMIT( FILES
    "${CTEST_SOURCE_DIRECTORY}/cmake/python/data/CDashSubprojectDependencies.xml"
    RETURN_VALUE SUBMIT_RETURN_VAL
    )
  MESSAGE("Submitted subproject dependencies: Return='${SUBMIT_RETURN_VAL}'")
  
  #
  # loop over all Trilinos packages
  #
  
  MESSAGE("\nBegin incremental building and testing of Trilinos pacakges ...\n")
  
  SET(Trilinos_LAST_WORKING_PACKAGE)
  SET(Trilinos_FAILED_PACKAGES)
  
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
  
    SET_PROPERTY(GLOBAL PROPERTY SubProject ${PACKAGE})
    SET_PROPERTY(GLOBAL PROPERTY Label ${PACKAGE})
  
    MESSAGE("\nCurrent Trilinos package ='${PACKAGE}'\n")
  
    #
    # Configure the package and its dependent packages
    #
  
    MESSAGE("\nConfiguring PACKAGE='${PACKAGE}'\n")
  
    # Create CONFIGURE_OPTIONS for this PACKAGE
    SET( CONFIGURE_OPTIONS
      "-DCTEST_USE_LAUNCHERS:BOOL=${CTEST_USE_LAUNCHERS}"
      "-DTrilinos_ENABLE_${PACKAGE}:BOOL=ON"
      "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
      "-DTrilinos_ENABLE_TESTS:BOOL=ON")
    IF (DEFINED Trilinos_LAST_WORKING_PACKAGE)
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_ENABLE_${Trilinos_LAST_WORKING_PACKAGE}:BOOL=")
      SET(Trilinos_LAST_WORKING_PACKAGE)
    ENDIF()
    FOREACH(FAILED_PACKAGE ${Trilinos_FAILED_PACKAGES})
      LIST(APPEND CONFIGURE_OPTIONS
        "-DTrilinos_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
    ENDFOREACH()
    MESSAGE("CONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")
  
    # Do the configure
    CTEST_CONFIGURE(
      BUILD "${CTEST_BINARY_DIRECTORY}"
      OPTIONS "${CONFIGURE_OPTIONS}" # New option!
      RETURN_VALUE CONFIGURE_RETURN_VAL
      )
  
    # If the configure failed, add the package to the list
    # of failed packages
    IF (NOT "${CONFIGURE_RETURN_VAL}" EQUAL "0")
      MESSAGE("${PACKAGE} FAILED to configure")
      LIST(APPEND Trilinos_FAILED_PACKAGES ${PACKAGE})
    ELSE()
      # load target properties and test keywords
      CTEST_READ_CUSTOM_FILES(BUILD "${CTEST_BINARY_DIRECTORY}")
    ENDIF()
  
    # Submit configure results and the notes to the dashboard 
    CTEST_SUBMIT( PARTS configure notes )
  
    #
    # If configure passed then try the build.  Otherwise, move on to
    # to the next package.
    #
  
    if ("${CONFIGURE_RETURN_VAL}" EQUAL "0")
  
      # Start by trying to build just the librries for the current package
  
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
  
      CTEST_SUBMIT( PARTS build APPEND )
  
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
        CTEST_SUBMIT( PARTS build APPEND )  
  
        IF (CTEST_DO_TEST)
          # Run the tests that match the ${PACKAGE} name 
          MESSAGE("\nRunning test for package '${PACKAGE}' ...\n")
          CTEST_TEST(BUILD "${CTEST_BINARY_DIRECTORY}"
            INCLUDE "^${PACKAGE}_"
            APPEND
            )
          CTEST_SUBMIT(PARTS Test APPEND)
        ENDIF()
  
        IF (CTEST_DO_COVERAGE_TESTING)
          MESSAGE("\nRunning coverage for package '${PACKAGE}' ...\n")
          CTEST_COVERAGE(BUILD "${CTEST_BINARY_DIRECTORY}")
          CTEST_SUBMIT(PARTS Test APPEND)
        ENDIF() 
 
        IF (CTEST_DO_MEMORY_TESTING)
          MESSAGE("\nRunning memory testing for package '${PACKAGE}' ...\n")
          CTEST_MEMCHECK(BUILD "${CTEST_BINARY_DIRECTORY}")
          CTEST_SUBMIT(PARTS Test APPEND)
        ENDIF()
  
        # Remember this package so we can turn it off
        SET(Trilinos_LAST_WORKING_PACKAGE "${PACKAGE}")
  
      ELSE()
  
        MESSAGE("FAILED library build for package '${PACKAGE}'")
        LIST(APPEND Trilinos_FAILED_PACKAGES ${PACKAGE})
  
      ENDIF()
  
    ENDIF()
  
  ENDFOREACH(PACKAGE)
  
  IF(Trilinos_FAILED_PACKAGES)
    MESSAGE("\nFinal set of failed packages: '${Trilinos_FAILED_PACKAGES}'")
  ENDIF()
  
  MESSAGE("\nDone with the incremental building and testing of Trilinos pacakges!\n")

ENDMACRO()
