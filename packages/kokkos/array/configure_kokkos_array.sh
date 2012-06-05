#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
# Additional command-line arguments given to this script will be
# passed directly to CMake.
#
EXTRA_ARGS=$@

#
# Force CMake to re-evaluate build options.
#
rm -rf CMakeCache.txt

# You will have to edit the paths below.
cmake \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_COMPILER:FILEPATH=/path/to/c++/compilier \
  -D CMAKE_C_COMPILER:FILEPATH=/path/to/c/compilier \
  -D CMAKE_Fortran_COMPILER:FILEPATH=/path/to/fortran \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -D CMAKE_CXX_FLAGS:STRING="-Wall -O3" \
  -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="-Werror" \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_Kokkos:BOOL=ON \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
  -D Kokkos_ENABLE_NodeAPI:BOOL=OFF \
  -D Kokkos_ENABLE_LinAlg:BOOL=OFF \
  -D Kokkos_ENABLE_Array:BOOL=ON \
  -D Kokkos_ENABLE_DeviceCuda:BOOL=ON \
  -D TPL_ENABLE_CUDA:BOOL=ON \
  -D TPL_ENABLE_Pthread:BOOL=ON \
  -D TPL_ENABLE_HWLOC:BOOL=ON \
  -D HWLOC_INCLUDE_DIRS:FILEPATH=/path/to/hwloc/inlude/directory/ \
  -D HWLOC_LIBRARY_DIRS:FILEPATH=/path/to/hwloc/libraries/directory \
  -D BLAS_INCLUDE_DIRS:FILEPATH=/path/to/blas/include/directory/ \
  -D BLAS_LIBRARY_DIRS:FILEPATH=/path/to/blas/libraries/directory/ \
  -D LAPACK_INCLUDE_DIRS:FILEPATH=/path/to/lapack/include/directory/ \
  -D LAPACK_LIBRARY_DIRS:FILEPATH=/path/to/lapack/libraries/dirctory/ \
  -D CUDA_NVCC_FLAGS:STRING="-arch=sm_20" \
  $EXTRA_ARGS \
  ../Trilinos



