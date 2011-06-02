
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test2.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME CONTINUOUS_${COMM_TYPE}_OPT_DEV_SHARED)
#SET(CTEST_TEST_TIMEOUT 900)

#override the default number of processors to run on.
SET( CTEST_BUILD_FLAGS "-j12 -i" )
SET( CTEST_PARALLEL_LEVEL "12" )

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

#disabling Mesquite because of a build error when shared libs is turned on.
SET(EXTRA_EXCLUDE_PACKAGES Mesquite STK Claps)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_ENABLE_DEBUG:BOOL=ON"
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  "-DMPI_BASE_DIR:PATH=/home/trilinos"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/include"
  "-DTPL_ENABLE_Expat:BOOL=ON"
  "-DTPL_ENABLE_LAMMPS:BOOL=ON"
  "-DTPL_ENABLE_couple:BOOL=ON"
  "-DTPL_ENABLE_SPPARKS:BOOL=ON"
  "-DTPL_ENABLE_SuperLU:BOOL=ON"
  )

# 2009/11/26: rabartl: Do we really want to be pointing to Trilinos_DATA_DIR?
# Unless the CI sever is going to be updatting this every iteration this is
# likely to cause an inconsistency and a test failure.  For example, if
# someone adds a new test and then adds new test data, the CI server will only
# get the new test but not the new test data.  Therefore, I am removing this.
# # "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"


#
# Set the rest of the system-specific options and run the dashboard build/test
#

SET(CTEST_TEST_TYPE Continuous)

SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY ON)
SET(CTEST_ENABLE_MODIFIED_PACKAGES_ONLY OFF)
SET(done 0)

WHILE(NOT done)
  SET(START_TIME ${CTEST_ELAPSED_TIME})

  TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

  SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY OFF)
  SET(CTEST_ENABLE_MODIFIED_PACKAGES_ONLY ON)

  MESSAGE("Before CTEST_SLEEP: CTEST_ELAPSED_TIME='${CTEST_ELAPSED_TIME}'")

  # Loop at most once every 3 minutes (180 seconds)
  CTEST_SLEEP(${START_TIME} 180 ${CTEST_ELAPSED_TIME})

  # Stop after 14 hours:
  IF(${CTEST_ELAPSED_TIME} GREATER 50400)
    SET(done 1)
  ENDIF()

  MESSAGE("Bottom of continuous while loop: CTEST_ELAPSED_TIME='${CTEST_ELAPSED_TIME}'")
ENDWHILE()
