INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.dash2win64.msvc.cmake")

#
# Set the options specific to this build case
#
#SET(CTEST_DO_UPDATES FALSE)
SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_RELEASE_PS)
#SET(CTEST_TEST_TIMEOUT 900)
#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
