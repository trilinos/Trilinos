
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
SET(BUILD_DIR_NAME ZOLTAN_PLUS_USERS)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=1800"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_MPI:BOOL=ON"
  "-DCMAKE_C_FLAGS:STRING=-std=c99"
  "-DMPI_EXEC_MAX_NUMPROCS:STRING=11"
  "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
  "-DTrilinos_ENABLE_EXAMPLES:BOOL=ON"
  "-DTrilinos_ENABLE_TESTS:BOOL=ON"
  "-DTrilinos_ENABLE_Zoltan:BOOL=ON"
  "-DZoltan_ENABLE_ParMETIS:BOOL=ON"
  "-DZoltan_ENABLE_Scotch:BOOL=ON"
  "-DParMETIS_INCLUDE_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-DParMETIS_LIBRARY_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-DScotch_INCLUDE_DIRS:FILEPATH=/home/lriesen/system/scotch_5.1.10a-32/include"
  "-DScotch_LIBRARY_DIRS:FILEPATH=/home/lriesen/system/scotch_5.1.10a-32/lib"
  "-DTrilinos_ENABLE_Isorropia:BOOL=ON"
  "-DTrilinos_ENABLE_Epetra:BOOL=ON"
  "-DTrilinos_ENABLE_EpetraExt:BOOL=ON"
  "-DTrilinos_ENABLE_ML:BOOL=ON"
  "-DTrilinos_ENABLE_STK:BOOL=ON"
  "-DTrilinos_ENABLE_Didasko:BOOL=ON"
  "-DTrilinos_ENABLE_trilinoscouplings:BOOL=ON"
  "-DTPL_BLAS_LIBRARIES:STRING=-L/usr/local/lib -lcblas -lf77blas -latlas -lblas"
  "-DTPL_LAPACK_LIBRARIES:STRING=/usr/local/lib/liblapack.a"
  )

SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )
SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
SET( CTEST_BUILD_FLAGS "-j4 -i" )
SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
SET( CTEST_MEMORYCHECK_COMMAND /usr/local/bin/valgrind )

TRILINOS_CTEST_DRIVER()
