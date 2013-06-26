#
# Primary Stable options for serial builds with hybrid GCC 4.6.1 C/C++ 
# Fortran
#

# Must be including first in order to define TRILINOS_TOOLSET_BASE
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.6.1-base-options.cmake)

# Set up the hybrid compilers
SET(CMAKE_CXX_COMPILER      "${GCC_BASE_DIR}/bin/g++"       CACHE FILEPATH "")
SET(CMAKE_C_COMPILER        "${GCC_BASE_DIR}/bin/gcc"       CACHE FILEPATH "")
SET(CMAKE_Fortran_COMPILER  "${GCC_BASE_DIR}/bin/gfortran"  CACHE FILEPATH "")
