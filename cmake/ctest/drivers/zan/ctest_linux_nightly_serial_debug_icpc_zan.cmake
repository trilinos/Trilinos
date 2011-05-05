
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.zan.icpc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG_ICPC)
SET(CTEST_TEST_TYPE Experimental)
SET(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:/opt/intel/Compiler/11.1/064/lib/intel64:/opt/intel/Compiler/11.1/064/mkl/lib/em64t")
#SET(CTEST_TEST_TIMEOUT 900)

#SET(Trilinos_PACKAGES Teuchos RTOp GlobiPack Thyra OptiPack Stratimikos Phalanx Rythmos)
SET(EXTRA_EXCLUDE_PACKAGES Didasko Sundance Piro Rythmos TrilinosCouplings NOX STK Pamgen Thyra Tpetra Zoltan Stokhos FEApp)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_Fortran:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
