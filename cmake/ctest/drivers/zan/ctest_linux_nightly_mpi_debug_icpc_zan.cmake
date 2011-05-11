
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.zan.icpc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME MPI_DEBUG_ICPC)
SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_TEST_TIMEOUT 900)
SET(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:/opt/intel/Compiler/11.1/064/lib/intel64:/opt/intel/Compiler/11.1/064/mkl/lib/em64t")

#SET(Trilinos_PACKAGES Teuchos RTOp GlobiPack Thyra OptiPack Stratimikos Phalanx Rythmos)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_Fortran:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
