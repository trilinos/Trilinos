
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Build Zoltan and all Trilinos packages that have optional or
#   required dependencies on it.
#
# Zoltan will use 32-bit global IDs.  ParMetis and
#    Scotch are using 32-bit global IDs.

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME ZOLTAN_PLUS_USERS)
SET(Trilinos_PACKAGES Zoltan)
SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)
SET(CTEST_TEST_TIMEOUT "1800")
#SET(CTEST_TEST_TYPE Nightly)
SET(CTEST_TEST_TYPE Experimental)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_MPI:BOOL=ON"
  "-DMPI_EXEC_MAX_NUMPROCS:STRING=11"
  "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
  "-DTrilinos_ENABLE_EXAMPLES:BOOL=ON"
  "-DTrilinos_ENABLE_TESTS:BOOL=ON"
  "-DTrilinos_ENABLE_Zoltan:BOOL=ON"
  "-DZoltan_ENABLE_UINT_IDS:BOOL=ON"
  "-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "-DZoltan_ENABLE_ParMETIS:BOOL=ON"
  "-DZoltan_ENABLE_Scotch:BOOL=ON"
  "-DParMETIS_INCLUDE_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-DParMETIS_LIBRARY_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-DScotch_INCLUDE_DIRS:FILEPATH=/home/lriesen/local/system/scotch_5.1.10a-32/include"
  "-DScotch_LIBRARY_DIRS:FILEPATH=/home/lriesen/local/system/scotch_5.1.10a-32/lib"
  "-DLAPACK_LIBRARY_DIRS:FILEPATH=/usr/local/lib"
  "-DBLAS_LIBRARY_DIRS:FILEPATH=/usr/local/lib"
  "-DTPL_LAPACK_LIBRARIES:STRING=-llapack -lblas -lgfortran"
  )

SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )
SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
SET( CTEST_BUILD_FLAGS "-j4 -i" )
SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
SET( CTEST_MEMORYCHECK_COMMAND /usr/local/bin/valgrind )

TRILINOS_CTEST_DRIVER()
