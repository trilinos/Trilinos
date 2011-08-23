# Primary Stable options for serial builds with GCC 4.5.1 C/C++ and Intel 12.0.4 Fortran
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-base-options.cmake)
# SET(CMAKE_Fortran_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/gfortran" CACHE FILEPATH "")
SET(CMAKE_Fortran_COMPILER "${INTEL_COMPILER_BASE}/bin/intel64/ifort" CACHE FILEPATH "")
SET(CMAKE_CXX_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/g++" CACHE FILEPATH "")
SET(CMAKE_C_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/gcc" CACHE FILEPATH "")
