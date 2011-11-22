#
# checkin-test-fissile4.sh MPI_DEBUG_SS build
#

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE  ON  CACHE BOOL  "")
SET(TPL_ENABLE_MPI  ON  CACHE  BOOL  "")
SET(Trilinos_ENABLE_CHECKED_STL  ON  CACHE  BOOL  "")
SET(Trilinos_ENABLE_DEBUG_SYMBOLS  ON  CACHE  BOOL  "")
SET(Teuchos_ENABLE_DEFAULT_STACKTRACE  OFF  CACHE  BOOL  "")

# Disable boost libs since not build with checked-stl turned on
SET(STK_ENABLE_BoostLib  OFF  CACHE  BOOL  "")

# Include last so that above override these cache variables
INCLUDE(${${PROJECT_NAME}_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-mpi-ss-options.cmake)
INCLUDE(${${PROJECT_NAME}_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-debug-options.cmake)
