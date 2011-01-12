
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.s909348.gcc.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../TrilinosVersion.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME "SERIAL_RELEASE_${Trilinos_VERSION}")
SET(Trilinos_TRACK ${Trilinos_TESTING_TRACK})
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_BRANCH ${Trilinos_REPOSITORY_BRANCH})

# Exclude Sundance because of strange segfault (see bug 4382)
SET(EXTRA_EXCLUDE_PACKAGES Sundance)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS=/Users/bmpersc/lib/boost_1_38_0"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
