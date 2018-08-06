# This file can be read in using either:
#
#   -C <abs-path>/<file-name>.cmake
#
# or:
#
#   -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/<file-name>.cmake
#

# Handle this being passed in with -C option instead of
# <Project>_CONFIGURE_OPTIONS_FILE.
IF ("${PROJECT_NAME}" STREQUAL "")
  SET(PROJECT_NAME Trilinos)
  INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../tribits/core/utils/AssertDefined.cmake")
ENDIF()

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/SetUtils.cmake")

TRIL_SET_BOOL_CACHE_VAR(${PROJECT_NAME}_ENABLE_OpenMP ON)

TRIL_SET_CACHE_VAR(MPI_EXEC_PRE_NUMPROCS_FLAGS "--bind-to;none"
  CACHE STSRING)
# NOTE: The above is a workaround for the problem of having threads on MPI
# ranks bind to the same cores (see #2422).

# Disable just one Teko sub-unit test that fails with GCC 4.8.4 + OpenMP (#2712)
TRIL_SET_BOOL_CACHE_VAR(Teko_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D TRUE)

# Disable these tests until they can get fixed (#2691)
TRIL_SET_BOOL_CACHE_VAR(ShyLU_DDFROSch_test_frosch_laplacian_epetra_2d_gdsw_MPI_4_DISABLE TRUE)
TRIL_SET_BOOL_CACHE_VAR(ShyLU_DDFROSch_test_frosch_laplacian_epetra_2d_rgdsw_MPI_4_DISABLE TRUE)
TRIL_SET_BOOL_CACHE_VAR(ShyLU_DDFROSch_test_frosch_interfacesets_2D_MPI_4_DISABLE TRUE)

# NOTE: The order of these includes matters!

include("${CMAKE_CURRENT_LIST_DIR}/MpiReleaseDebugSharedPtSettings.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/BasicCiTestingSettings.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/sems/SEMSDevEnv.cmake")
