INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.kallikrates.msvc.cmake")

#
# Set the options specific to this build case
#
#SET(CTEST_DO_UPDATES FALSE)
SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_RELEASE)
#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
