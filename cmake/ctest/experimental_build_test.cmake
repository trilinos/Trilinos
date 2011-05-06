#
# CTest script that is used to do an experimental build/test right from a
# developer's own build directory.
#
# NOTE: You need a recent (CVS) version of CMake/CTest for this to work.  You
# can't use CMake/CTest 2.6.x.  If you don't have the right version, it will
# tell you.
#
# To run this script:
#
# 1) First configure your directory without enabling any packages (to set the
# platform-specific and site-specific options).  The CMakeCache.txt file that
# is created will then be modified by the script as packages are enabled and
# disabled.  NOTE: For the safest results, start with an empty build directory
# (except for your configure script of course).  This seems to be needed when
# switching back and forth between an 'Experimental' and 'Nightly' build for
# instance.
#
# 2) Run the script (overriding any appropriate options) as:
#
#    env Trilinos_PACKAGES="<PACKAGES>" \
#      ctest -S $TRILINOS_HOME/cmake/ctest/experimental_build_test.cmake -VV
#
# where PACAKGES is the semi-colon-separated list of packages being tested
# (e.g. Trilinos_PACKAGES="Teuchos;Epetra;NOX") and TRILINOS_HOME points back
# to your home Trilinos directory.  You can take off the -VV argument if you
# don't want this to be too verbose.
#
# There are a number of other options that you can change as
# environment varibles.  See the macros SET_DEFAULT_AND_FROM_ENV(...)
# in the file TrilinosCTestDriverCore.cmake.  One option that you
# might want to overridde, for instance is CTEST_BUILD_NAME so that
# you can insert a special name into the dashboard.
#
# When this script finishes running, the last package listed in
# Trilinos_PACAKGES will be enabled in the CMakeCache.txt file.
#
# NOTE: It is better to use the CMake-built make target 'experimental' to run
# this script as it takes care of the details of manipulating the cache and
# restoring the package enables when it is done.
#


#
# General setup code:
#
# Do not modify any of this directly, use use environment variables instead!
#

#
# Include some CMake/CTest code files
#

SET( CMAKE_MODULE_PATH
  "${CTEST_SCRIPT_DIRECTORY}"
  "${CTEST_SCRIPT_DIRECTORY}/../utils"
  )

INCLUDE(TrilinosCTestDriverCore)
INCLUDE(GetLastDirName)

#
# Override some configuration variables
#

# All these can be changed by env vars
SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_DO_UPDATES FALSE)
SET(Trilinos_WARNINGS_AS_ERRORS_FLAGS "-Werror")

# Don't change these in the env!
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)
SET(CTEST_GENERATE_DEPS_XML_OUTPUT_FILE TRUE)
SET(CTEST_WIPE_CACHE FALSE)

SET(CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}/../..")
GET_FILENAME_COMPONENT(PWD . REALPATH)
SET(CTEST_BINARY_DIRECTORY "${PWD}")
SET(CTEST_NOTES_FILES "${CTEST_BINARY_DIRECTORY}/do-configure")

GET_LAST_DIR_NAME("${CTEST_BINARY_DIRECTORY}" BUILD_DIR_NAME)

# Can be overridden by the environment
SET( CTEST_BUILD_NAME "${HOST_TYPE}-${BUILD_DIR_NAME}" )
SET( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES OFF )

#
# Run the build/test/submit driver
#

TRILINOS_CTEST_DRIVER()
