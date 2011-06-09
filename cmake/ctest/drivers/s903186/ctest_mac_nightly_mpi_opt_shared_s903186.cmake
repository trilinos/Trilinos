
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.s903186.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV_SHARED)
#SET(CTEST_TEST_TIMEOUT 900)

SET(EXTRA_EXCLUDE_PACKAGES Zoltan STK Mesquite Claps FEI Phdmesh NOX Piro Sundance)

#Disabling STK below also so it is not turned on by TrilinosCouplings
SET( EXTRA_CONFIGURE_OPTIONS
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  "-DSWIG_EXECUTABLE:FILEPATH=/Users/jmwille/install/bin/swig"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DTrilinos_ENABLE_STK:BOOL=OFF"
  "-DNOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
