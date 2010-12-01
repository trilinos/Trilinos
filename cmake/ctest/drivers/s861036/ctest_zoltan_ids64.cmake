
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#   Build Zoltan, Isorropia, Epetra and EpetraExt
#   Specify 64-bit global IDs for Zoltan and Isorropia
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=120"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-D CMAKE_C_FLAGS:STRING=-std=c99 -DZOLTAN_ID_TYPE_LONG"
  "-D CMAKE_CXX_FLAGS:STRING=-DZOLTAN_ID_TYPE_LONG"
  "-D MPI_EXEC_MAX_NUMPROCS:STRING=11"
  "-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
  "-D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
  "-D Trilinos_ENABLE_TESTS:BOOL=ON"
  "-D Trilinos_ENABLE_Zoltan:BOOL=ON"
  "-D Trilinos_ENABLE_Isorropia:BOOL=ON"
  "-D Trilinos_ENABLE_Epetra:BOOL=ON"
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON"
  "-D Zoltan_ENABLE_ParMETIS:BOOL=ON"
  "-D Zoltan_ENABLE_Scotch:BOOL=ON"
  "-D ParMETIS_INCLUDE_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-D ParMETIS_LIBRARY_DIRS:FILEPATH=/home/lafisk/system/parmetis/ParMetis-3.1"
  "-D Scotch_INCLUDE_DIRS:FILEPATH=/home/lriesen/system/scotch_5.1.10a-32/include"
  "-D Scotch_LIBRARY_DIRS:FILEPATH=/home/lriesen/system/scotch_5.1.10a-32/lib"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
