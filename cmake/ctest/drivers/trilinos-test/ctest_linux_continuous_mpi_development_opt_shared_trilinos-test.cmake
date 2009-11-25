
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME CONTINUOUS_${COMM_TYPE}_OPT_DEV_SHARED)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

#disabling Mesquite because of a build error when shared libs is turned on.
SET(EXTRA_EXCLUDE_PACKAGES Mesquite)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=120"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  "-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost_1_38_0"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

SET(CTEST_TEST_TYPE Continuous)

SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY ON)
SET(done 0)

WHILE(NOT done)
  SET(START_TIME ${CTEST_ELAPSED_TIME})

  TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

  SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY OFF)

  MESSAGE("Before CTEST_SLEEP: CTEST_ELAPSED_TIME='${CTEST_ELAPSED_TIME}'")

  # Loop at most once every 3 minutes (180 seconds)
  CTEST_SLEEP(${START_TIME} 180 ${CTEST_ELAPSED_TIME})

  # Stop after 10 hours:
  IF(${CTEST_ELAPSED_TIME} GREATER 36000)
    SET(done 1)
  ENDIF()

  MESSAGE("Bottom of continuous while loop: CTEST_ELAPSED_TIME='${CTEST_ELAPSED_TIME}'")
ENDWHILE()
