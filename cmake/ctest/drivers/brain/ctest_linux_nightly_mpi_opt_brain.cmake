
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../TrilinosVersion.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.brain.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME "MPI_OPT_DEV")
SET(Trilinos_TRACK ${Trilinos_TESTING_TRACK})
SET(Trilinos_BRANCH ${Trilinos_REPOSITORY_BRANCH})

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

