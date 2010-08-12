

#
# General setup code:
#
# Do not modify any of this directly, use use environment variables instead!
#

#
# Include some CMake/CTest code files
#

#MESSAGE("CTEST_SCRIPT_DIRECTORY = ${CTEST_SCRIPT_DIRECTORY}")

SET( CMAKE_MODULE_PATH
  "${CTEST_SCRIPT_DIRECTORY}/.."
  "${CTEST_SCRIPT_DIRECTORY}/../../utils"
  )


#
# Includes
#

INCLUDE(TrilinosCTestDriverCore)
INCLUDE(GetLastDirName)


#
# Override some configuration variables
#

# All these can be changed by env vars
SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_DO_UPDATES FALSE)
SET(Trilinos_WARNINGS_AS_ERRORS_FLAGS "-Werror")

# Don't change these in the env!
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)
SET(CTEST_GENERATE_DEPS_XML_OUTPUT_FILE TRUE)
SET(CTEST_WIPE_CACHE FALSE)

SET(CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}/../../..")
SET(Trilinos_DEPS_HOME_DIR "${CTEST_SOURCE_DIRECTORY}/cmake/DependencyUnitTests/MockTrilinos")
GET_FILENAME_COMPONENT(PWD . REALPATH)
SET(CTEST_BINARY_DIRECTORY "${PWD}")

SET_DEFAULT_AND_FROM_ENV(Trilinos_CTEST_COMMAND ctest)
SET(CTEST_COMMAND ${Trilinos_CTEST_COMMAND})


#
# Run the build/test/submit driver
#

TRILINOS_CTEST_DRIVER()
