#-------------------------------------------------------------------------------
# Nightly testing on linux platform thumper.mp.sandia.gov
# Debug with Coverage and MemoryCheck
#-------------------------------------------------------------------------------

SET(CTEST_SOURCE_NAME Trilinos)
SET(TEST_TYPE nightly)
SET(BUILD_TYPE debug)
SET(EXTRA_BUILD_TYPE serial)
SET(HOSTTYPE Linux) # Have to set this manually on this machine for some reason?

SET(CTEST_DASHBOARD_ROOT /home/rabartl/PROJECTS/dashboards/Trilinos/PACKAGE_DEPS)
SET(CTEST_CMAKE_COMMAND /usr/local/bin/cmake)

# Options for Nightly builds
#SET(CTEST_BACKUP_AND_RESTORE TRUE)
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)
#SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)
SET(CTEST_CVS_CHECKOUT
  "cvs -q -d :ext:software.sandia.gov:/space/CVS co ${CTEST_SOURCE_NAME}"
)
SET (CTEST_CVS_COMMAND
  "cvs -q -d :ext:software.sandia.gov:/space/CVS co ${CTEST_SOURCE_NAME}"
)

SET(CTEST_BINARY_NAME BUILD)
SET(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
SET(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

SET(CTEST_COMMAND 
  "\"${CTEST_EXECUTABLE_NAME}\" -D NightlyStart"
  "\"${CTEST_EXECUTABLE_NAME}\" -D NightlyConfigure"
  "\"${CTEST_EXECUTABLE_NAME}\" -D NightlyBuild"
  "\"${CTEST_EXECUTABLE_NAME}\" -D NightlySubmit"
  "\"${CTEST_EXECUTABLE_NAME}\" -D NightlyTest"
  "\"${CTEST_EXECUTABLE_NAME}\" -D NightlySubmit -A \"${CTEST_BINARY_DIRECTORY}/CMakeCache.txt\""
)

SET(CTEST_INITIAL_CACHE "

Trilinos_ENABLE_CXX:BOOL=OFF
CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE
Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=ON
Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF

BUILDNAME:STRING=${HOSTTYPE}-TrilinosDepUnitTests

CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}

MAKECOMMAND:STRING=gmake -j 8

")
