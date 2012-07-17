#
# Primary Stable options for MPI builds with hybrid GCC 4.6.1 C/C++ 
# Fortran
#

# Included first ot define TRILINOS_TOOLSET_BASE and MKLROOT
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.6.1-base-options.cmake)

# Point to OpenMPI build with GCC 4.6.1 C/C++ 
SET(MPI_BASE_DIR "${TRILINOS_TOOLSET_BASE}" CACHE PATH "")
