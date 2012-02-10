#
# Primary Stable options for MPI builds with hybrid GCC 4.6.1 C/C++ 
# Fortran
#

# Included first ot define TRILINOS_TOOLSET_BASE and MKLROOT
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.6.1-base-options.cmake)

# Point to OpenMPI build with GCC 4.6.1 C/C++ 
SET(MPI_BASE_DIR "${TRILINOS_TOOLSET_BASE}" CACHE PATH "")

# Used only when TriKota/Dakota are enabled with CMake
SET(TriKota_ENABLE_DakotaCMake ON CACHE BOOL "")
SET(ENABLE_DAKOTA_TESTS OFF CACHE BOOL "")
SET(HAVE_ACRO OFF CACHE BOOL "")
SET(HAVE_AMPL OFF CACHE BOOL "")
SET(HAVE_X_GRAPHICS OFF CACHE BOOL "")
SET(HAVE_HOPSPACK OFF CACHE BOOL "")
