INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.kallikrates.msvc.cmake")

#
# Set the options specific to this build case
#
#SET(CTEST_DO_UPDATES FALSE)
SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_OPT_DEV_SHARED)
SET(CTEST_TEST_TYPE EXPERIMENTAL)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

#The dll export macros have only been added for a few  packages so we can really only test those packages
SET(Trilinos_PACKAGES "Teuchos;Epetra;Anasazi")

SET( EXTRA_CONFIGURE_OPTIONS
  "-DBUILD_SHARED_LIBS:BOOL=ON"
)

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
