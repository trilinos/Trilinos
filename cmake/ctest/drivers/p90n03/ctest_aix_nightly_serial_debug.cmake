INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.p90n03.xl.cmake")

#
# Set the options specific to this build case
#
#SET(CTEST_DO_UPDATES FALSE)
SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG)
#SET(CTEST_TEST_TIMEOUT 900)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
