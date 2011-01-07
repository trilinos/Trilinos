
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.gabriel.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG_NOTRACE)
#SET(CTEST_TEST_TIMEOUT 900)

# Only test packages that uses RCP up to Thyra
SET(Trilinos_PACKAGES Teuchos RTOp EpetraExt Thyra)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON"
  "-DTeuchos_ENABLE_DEBUG_RCP_NODE_TRACING:BOOL=OFF"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
