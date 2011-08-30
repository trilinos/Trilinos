# Primary Stable options for serial builds with GCC 4.6.1
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.6.1-base-options.cmake)
SET(CMAKE_Fortran_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/gfortran" CACHE FILEPATH "")
SET(CMAKE_CXX_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/g++" CACHE FILEPATH "")
SET(CMAKE_C_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/gcc" CACHE FILEPATH "")
