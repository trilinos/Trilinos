
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.u233.icpc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_RELEASE_ICPC)
#SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_PACKAGES Teuchos RTOp GlobiPack Tpetra Thyra OptiPack Stratimikos Rythmos MOOCHO)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTPL_ENABLE_Boost:BOOL=OFF"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
