
#
# Set the options specific to this build
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME PACKAGE_DEPS)


SET(CTEST_INITIAL_CACHE "

CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE

Trilinos_ENABLE_CXX:BOOL=OFF
Trilinos_ENABLE_C:BOOL=OFF
Trilinos_ENABLE_Fortran:BOOL=OFF

Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=ON

")


#
# Read in the platform-independent and platform-dependent options
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestSupport.godel.cmake")
