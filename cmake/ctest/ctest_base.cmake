#
# Default varaible settings
#

SET(CTEST_SOURCE_NAME Trilinos)
SET(CTEST_BINARY_NAME BUILD)

SET( CTEST_DASHBOARD_ROOT
  "${CTEST_SCRIPT_DIRECTORY}/../${BUILD_DIR_NAME}"
  )

SET( CTEST_SOURCE_DIRECTORY
  "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}"
  )
SET( CTEST_BINARY_DIRECTORY
  "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

SET( CTEST_NOTES_FILES_DEFAULT
  "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt"
  "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}"
  )

SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")

SET(DEFAULT_TEST_TYPE Experimental)

SET(DEFAULT_CTEST_DROP_METHOD http)

FIND_PROGRAM(HOSTNAME_EXE NAMES hostname)
EXECUTE_PROCESS(
  COMMAND ${HOSTNAME_EXE}
  OUTPUT_VARIABLE CTEST_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

SET(CTEST_DROP_SITE "trilinos.sandia.gov")
SET(DEFAULT_CTEST_DROP_LOCATION "/cdash/submit.php?project=Trilinos")
SET(DEFAULT_CTEST_TRIGGER_SITE "")

SET(CTEST_ENABLE_COVERAGE OFF)
SET(CTEST_ENABLE_MEMCHECK OFF)


#
# Macro to run after all variables are set
#

# ToDo:
#
# (*) Put in support for optional coverage testing
#
# (*) Put in support for optional memory testing

MACRO(DO_TRILINOS_CTEST)

  SET(ENV_TEST_TYPE $ENV{CTEST_TEST_TYPE})
  IF (ENV_TEST_TYPE)
    SET(TEST_TYPE ${ENV_TEST_TYPE})
  ELSE()
    SET(TEST_TYPE ${DEFAULT_TEST_TYPE})
  ENDIF()

  SET(ENV_LEAVE_BINARY_DIR $ENV{CTEST_LEAVE_BINARY_DIR})

  IF (NOT ENV_LEAVE_BINARY_DIR)
    MESSAGE("\nEmptying the binary directory ...\n")
    CTEST_EMPTY_BINARY_DIRECTORY("${CTEST_BINARY_DIRECTORY}")
  ENDIF()

  MESSAGE("\nStarting the dashboard entry ...\n")
  CTEST_START(${TEST_TYPE})

  MESSAGE("\nUpdating the sources ...\n")
  CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}")

  MESSAGE("\nConfiguring the binary directory with cmake ...\n")
  FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt"
"${CTEST_INITIAL_CACHE}
DROP_METHOD:STRING=${DEFAULT_CTEST_DROP_METHOD}"
   )
  CTEST_CONFIGURE(BUILD "${CTEST_BINARY_DIRECTORY}")

  MESSAGE("\nReading custom files ...\n")
  CTEST_READ_CUSTOM_FILES("${CTEST_BINARY_DIRECTORY}")

  MESSAGE("\nBuilding ...\n")
  CTEST_BUILD(BUILD "${CTEST_BINARY_DIRECTORY}")

  MESSAGE("\nSumitting update, configure, and build results ...\n")
  CTEST_SUBMIT()
  # 20008/12/18: rabartl: Above, it looks like the CDash dashboard
  # only supports submitting after the build because it puts in
  # an entry after the build anyway and it will show a build with
  # no errors and no warnings if you submit before the build.

  MESSAGE("\nRunning the tests ...\n")
  CTEST_TEST(BUILD "${CTEST_BINARY_DIRECTORY}")

  MESSAGE("\nSubmitting test results and the notes files ...\n")
  SET(CTEST_NOTES_FILES ${CTEST_NOTES_FILES_DEFAULT})
  CTEST_SUBMIT()

  IF (CTEST_ENABLE_COVERAGE)

    MESSAGE("\nPerforming coverage testing ...\n")
    CTEST_COVERAGE(BUILD "${CTEST_BINARY_DIRECTORY}")

    MESSAGE("\nSumitting coverage results ...\n")
    CTEST_SUBMIT()

  ENDIF()

  IF (CTEST_ENABLE_MEMCHECK)

    MESSAGE("\nPerforming memory testing ...\n")
    CTEST_MEMCHECK(BUILD "${CTEST_BINARY_DIRECTORY}")

    MESSAGE("\nSumitting memory testing results ...\n")
    CTEST_SUBMIT()

  ENDIF()

ENDMACRO()
