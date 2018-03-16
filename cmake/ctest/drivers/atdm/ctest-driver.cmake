################################################################################
#
# Common CTest -S script for all ATDM builds
#
################################################################################

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.atdm.cmake")

# Run the genetic ATDM driver
TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
