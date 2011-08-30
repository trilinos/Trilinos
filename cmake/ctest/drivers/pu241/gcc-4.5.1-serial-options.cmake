# Primary Stable options for serial builds with GCC 4.5.1 C/C++ and Intel 12.0.4 Fortran
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-base-options.cmake)
# SET(CMAKE_Fortran_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/gfortran" CACHE FILEPATH "")
# SET(BLAS_INCLUDE_DIRS   ${MKLROOT}/include/gcc-4.5.1                     CACHE PATH     "Path to MKL BLAS Fortran modules compatible with GCC 4.5.1")
# SET(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS}                             CACHE PATH     "Path to MKL LAPACK Fortran modules compatible with GCC 4.5.1")

SET(CMAKE_Fortran_COMPILER "${INTEL_COMPILER_BASE}/bin/intel64/ifort" CACHE FILEPATH "")
SET(BLAS_INCLUDE_DIRS   ${MKLROOT}/include/intel64/lp64                     CACHE PATH  "Path to MKL BLAS Fortran modules compatible with Intel fortran")
SET(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS}                                CACHE PATH  "Path to MKL LAPACK Fortran modules compatible with Intel fortran")

SET(CMAKE_CXX_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/g++" CACHE FILEPATH "")
SET(CMAKE_C_COMPILER "${TRILINOS_TOOLSET_BASE}/bin/gcc" CACHE FILEPATH "")
