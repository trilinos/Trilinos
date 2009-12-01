
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME MPI_DEBUG_CONTINUOUS)
SET(CTEST_TEST_TYPE CONTINUOUS)
#SET(CTEST_TEST_TYPE EXPERIMENTAL)
# Wipe binary tree once per script invocation
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE TRUE)
# How long to run in minutes
SET (CTEST_CONTINUOUS_DURATION 600)
SET (CTEST_CONTINUOUS_INTERVAL 20)
SET (CTEST_TESTING_TIMEOUT 720)
#SET(CTEST_DO_COVERAGE_TESTING TRUE)
#SET(CTEST_DO_MEMORY_TESTING TRUE)
#override the default number of processors to run on.
SET( CTEST_BUILD_FLAGS "-j11 -i" )
SET( CTEST_PARALLEL_LEVEL "11" )

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCTEST_TESTING_TIMEOUT:STRING=900"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  )
#  "-DTPL_ENABLE_ParMETIS:BOOL=ON"
#  "-DParMETIS_LIBRARY_DIRS:PATH=/home/kddevin/code/ParMETIS3_1"
#  "-DTPL_ENABLE_Scotch:BOOL=ON"
#  "-DScotch_INCLUDE_DIRS:PATH=/home/kddevin/code/scotch_5.1/include"
#  "-DScotch_LIBRARY_DIRS:PATH=/home/kddevin/code/scotch_5.1/lib"
#  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
