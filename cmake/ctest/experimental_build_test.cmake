#
# CTest script that is used to do an experimental build/test right from a
# developer's own build directory.
#
# To run this script:
#
# 1) First configure your directory without enabling any packages (to
# set the platform-specific and site-specific options).  The
# CMakeCache.txt file that is created will then be modified by the
# script as packages are enabled and disabled.
#
# 2) Set the environment varaible TRILINOS_HOME to the base Trilinos source
# directory.  Or, you can set it in the 'env' command below.
#
# 3) Run the script (overriding any appropriate options) as:
#
#    env \
#        TRILINOS_HOME=../../../Trilinos \
#        Trilinos_PACKAGES="<PACKAGES>" \
#      ctest -S $TRILINOS_HOME/cmake/ctest/experimental_build_test.cmake \
#        -VV
#
# where PACAKGES is the semi-colon-separated list of packages being
# tested (e.g. Trilinos_PACKAGES="Teuchos;Epetra;NOX").  You can take
# off the -VV argument if you don't want this to be so verbose.
#
# There are a number of other options that you can change as
# environment varibles.  See the macros SET_DEFAULT_AND_FROM_ENV(...)
# in the file TrilinosCTestDriverCore.cmake.  One option that you
# might want to overridde, for instance is CTEST_BUILD_NAME so that
# you can insert a special name into the dashboard.
#
# When this script finishes running, the last package will be enabled
# in the CMakeCache.txt file.
#

SET( CMAKE_MODULE_PATH
  "${CTEST_SCRIPT_DIRECTORY}"
  "${CTEST_SCRIPT_DIRECTORY}/../utils"
  )

INCLUDE(TrilinosCTestDriverCore)
INCLUDE(GetLastDirName)

SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_WIPE_CACHE FALSE)
SET(CTEST_DO_UPDATES FALSE)
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

SET(CTEST_SOURCE_DIRECTORY "$ENV{TRILINOS_HOME}")
SET(CTEST_BINARY_DIRECTORY "$ENV{PWD}")

SET(CTEST_NOTES_FILES "${CTEST_BINARY_DIRECTORY}/do-configure")

GET_LAST_DIR_NAME("${CTEST_BINARY_DIRECTORY}"  BUILD_DIR_NAME)

SET_DEFAULT_AND_FROM_ENV( CTEST_BUILD_NAME "${HOST_TYPE}-${BUILD_DIR_NAME}" )

TRILINOS_CTEST_DRIVER()
