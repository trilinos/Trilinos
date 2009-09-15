
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_RELEASE_10.0)
SET(Trilinos_TRACK "Nightly Release 10.0")

SET(Trilinos_BRANCH "-r trilinos-release-10-0-branch")

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
  "-DDART_TESTING_TIMEOUT:STRING=120"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
