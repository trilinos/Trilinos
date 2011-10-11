
# Secondary Stable MPI build with GCC 4.5.1

# Must turn off this option so that Panzer/Drekar can work with Intrepid
SET(Intrepid_ENABLE_DEBUG_INF_CHECK OFF CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-mpi-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/boost-1.46.1-options.cmake)
