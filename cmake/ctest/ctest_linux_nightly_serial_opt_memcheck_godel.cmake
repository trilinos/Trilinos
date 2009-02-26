
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.godel.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_RELEASE_MEMCHECK)

SET(CTEST_DO_TEST FALSE)
SET(CTEST_DO_MEMORY_TESTING TRUE)

SET(Trilinos_PACKAGES Teuchos RTOp Epetra GlobiPack Tpetra EpetraExt
  Sacado Thyra OptiPack AztecOO Ifpack ML Stratimikos Rythmos MOOCHO)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=600"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
