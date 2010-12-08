
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Build Zoltan and all Trilinos packages that have optional or
#   required dependencies on it.
#
# Zoltan will use 32-bit global IDs.  ParMetis and
#    Scotch are using 32-bit global IDs.
#
# Didsko, Epetra EpetraExt Isorropia ML STK TrilinosCouplings
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_ZOLTAN_PLUS)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=300"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_MPI:BOOL=ON"
  "-D CMAKE_C_FLAGS:STRING=-std=c99"
  "-D MPI_EXEC_MAX_NUMPROCS:STRING=11"
  "-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
  "-D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
  "-D Trilinos_ENABLE_TESTS:BOOL=ON"
  "-D Trilinos_ENABLE_Zoltan:BOOL=ON"
  "-D Zoltan_ENABLE_ParMETIS:BOOL=ON"
  "-D Zoltan_ENABLE_Scotch:BOOL=ON"
  "-D ParMETIS_INCLUDE_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-D ParMETIS_LIBRARY_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-D Scotch_INCLUDE_DIRS:FILEPATH=/home/lriesen/system/scotch_5.1.10a-32/include"
  "-D Scotch_LIBRARY_DIRS:FILEPATH=/home/lriesen/system/scotch_5.1.10a-32/lib"
  "-D Trilinos_ENABLE_Isorropia:BOOL=ON"
  "-D Trilinos_ENABLE_Epetra:BOOL=ON"
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON"
  "-D Trilinos_ENABLE_ML:BOOL=ON"
  "-D Trilinos_ENABLE_STK:BOOL=ON"
  "-D Trilinos_ENABLE_Didasko:BOOL=ON"
  "-D Trilinos_ENABLE_trilinoscouplings:BOOL=ON"
  )

SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )
SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
SET( CTEST_BUILD_FLAGS "-j8 -i" )
SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
SET( CTEST_MEMORYCHECK_COMMAND /usr/local/bin/valgrind )

TRILINOS_CTEST_DRIVER()
