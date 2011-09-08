
# Include the base Windows configuration.
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/Windows.cmake)

# Additional parameters for the MPI build.
set(TPL_ENABLE_MPI ON CACHE BOOL "")
set(MPI_BASE_DIR "C:/Program Files (x86)/MPICH2")
set( Trilinos_ENABLE_Mesquite OFF CACHE BOOL "" )

