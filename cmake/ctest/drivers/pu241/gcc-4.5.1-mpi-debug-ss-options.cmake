# Secondary Stable MPI DEBUG build with GCC 4.5.1 (same as MPI_DEBUG_SS in checkin-test-fissile4.sh)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-mpi-ss-options.cmake)
SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE  ON  CACHE BOOL  "")
SET(CMAKE_BUILD_TYPE  DEBUG  CACHE  STRING  "")
SET(TPL_ENABLE_MPI  ON  CACHE BOOL  "")
SET(Trilinos_ENABLE_CHECKED_STL  ON  CACHE BOOL  "")
SET(Trilinos_ENABLE_DEBUG_SYMBOLS  ON  CACHE  BOOL  "")
SET(euchos_ENABLE_DEFAULT_STACKTRACE  OFF  CACHE  BOOL  "")

# ToDo: Re-enable stack tracing for nightly builds?
