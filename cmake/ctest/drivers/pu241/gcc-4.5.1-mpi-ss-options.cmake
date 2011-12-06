
# Secondary Stable MPI build with GCC 4.5.1

# Must turn off this option so that Panzer/Drekar can work with Intrepid
SET(Intrepid_ENABLE_DEBUG_INF_CHECK OFF CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.5.1-mpi-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/boost-1.46.1-options.cmake)
