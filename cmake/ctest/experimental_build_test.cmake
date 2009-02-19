#
# CTest script that is used to do an experimental build/test right from a
# developer's own build directory.
#
# To run this script:
#
# 1) First configure your directory without enabling any packages (to set
# the site-specific options).  The CMakeCache.txt file created will then be
# modified by the script as packages are enabled and disabled.
#
# 2) Set the environment varaible TRILINOS_HOME to the base Trilinos source
# directory.  Or, you can set it in the 'env' command below
#
# 3) Run the script (overriding any appropriate options) as:
#
#    env  Trilinos_PACKAGES="<PACKAGES>" \
#      ctest -S $TRILINOS_HOME/cmake/ctest/experimental_build_test.cmake \
#      -VV
#
# where PACAKGES is the semi-colon-separated list of packages being tested.
# You can take off the -VV argument if you don't want this to be so verbose.
#

INCLUDE(${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.cmake)

SET(CTEST_TEST_TYPE Experimental)

SET(CTEST_SOURCE_DIRECTORY "$ENV{TRILINOS_HOME}")
SET(CTEST_BINARY_DIRECTORY "$ENV{PWD}")

SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")

SET(CTEST_NOTES_FILES "${CTEST_BINARY_DIRECTORY}/do-configure")

# Get the name of the build directory we are in
FILE(TO_CMAKE_PATH ${CTEST_BINARY_DIRECTORY} STANDARD_CTEST_BINARY_DIRECTORY)
STRING(REGEX REPLACE "/.+/(.+)" "\\1" BUILD_DIR_NAME "${STANDARD_CTEST_BINARY_DIRECTORY}")

SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_NAME "${HOST_TYPE}-${BUILD_DIR_NAME}" )

SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

TRILINOS_CTEST_DRIVER()
