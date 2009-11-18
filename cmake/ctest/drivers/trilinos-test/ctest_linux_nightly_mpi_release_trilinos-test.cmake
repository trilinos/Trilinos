
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_RELEASE_10.0)
SET(Trilinos_TRACK "Nightly Release 10.0")

SET(Trilinos_BRANCH "trilinos-release-10-0-branch")

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=120"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost_1_38_0"
  "-DTPL_ENABLE_ParMETIS:BOOL=ON"
  "-DParMETIS_LIBRARY_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/ParMETIS3_1"
  "-DTPL_ENABLE_Scotch:BOOL=ON"
  "-DScotch_INCLUDE_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/scotch_5.1/include"
  "-DScotch_LIBRARY_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/scotch_5.1/lib "
  "-DTPL_ENABLE_ExodusII:BOOL=ON"
  "-DTPL_ENABLE_Nemesis:BOOL=ON"
  "-DNemesis_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/nemesis_3.09/include"
  "-DNemesis_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/nemesis_3.09/lib"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
