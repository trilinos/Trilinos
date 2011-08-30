
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test2.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV)
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DMPI_BASE_DIR:PATH=/home/trilinos"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/include"
  "-DTPL_ENABLE_ParMETIS:BOOL=ON"
  "-DParMETIS_LIBRARY_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/openmpi1.4/ParMETIS_3.1"
  "-DTPL_ENABLE_Scotch:BOOL=ON"
  "-DScotch_INCLUDE_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/openmpi1.4/Scotch_5.1.7/include"
  "-DScotch_LIBRARY_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/openmpi1.4/Scotch_5.1.7/lib "
  "-DTPL_ENABLE_HDF5:BOOL=ON"
  "-DHDF5_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/include"
  "-DHDF5_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/lib"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
