#
# Builds TriKota for a Trilinos bound build
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.pu241.gcc.4.5.1.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/SubmitToTrilinos.cmake")

SET(Trilinos_EXTRAREPOS_FILE
  "${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/ExtraExternalRepositories.dakota.cmake")

SET(COMM_TYPE MPI)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME MPI_DEBUG_GCC_TRIKOTA)
#SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_TEST_TIMEOUT 900)
SET(Trilinos_PACKAGES TriKota Piro LIME)
SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
SET(EXTRA_CONFIGURE_OPTIONS
  -DTrilinos_ENABLE_TriKota:BOOL=
  )

# Above, we override -DTrilinos_ENABLE_TriKota:BOOL=OFF in the default
# extra configure options but now -DTrilinos_ENABLE_TriKota:BOOL=ON
# will be set when processing the package TriKota

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
