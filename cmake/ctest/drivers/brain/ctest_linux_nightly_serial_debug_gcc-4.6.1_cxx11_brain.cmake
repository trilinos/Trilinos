
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../TrilinosVersion.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.brain.gcc-4.6.1.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME "SERIAL_DEBUG_GCC-4.6.1_CXX11_DEV")
SET(COMPILER_VERSION "GCC-4.6.1-C++1X")
SET(Trilinos_TRACK ${Trilinos_TESTING_TRACK})
SET(Trilinos_BRANCH ${Trilinos_REPOSITORY_BRANCH})

SET( EXTRA_CONFIGURE_OPTIONS
  "${EXTRA_CONFIGURE_OPTIONS}"
  "-DTrilinos_ENABLE_CXX11:BOOL=ON"
  "-DCMAKE_CXX_FLAGS:STRING=-std=c++0x"
)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

