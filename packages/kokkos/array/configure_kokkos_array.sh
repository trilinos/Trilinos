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
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -D CMAKE_CXX_FLAGS:STRING="-Wall -O3" \
  -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="-Werror" \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_KokkosArray:BOOL=ON \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
  -D KokkosArray_ENABLE_DeviceCuda:BOOL=ON \
  -D TPL_ENABLE_CUDA:BOOL=ON \
  -D TPL_ENABLE_CUSPARSE:BOOL=ON \
  -D CUDA_VERBOSE_BUILD:BOOL=ON \
  -D TPL_ENABLE_Pthread:BOOL=ON \
  -D TPL_ENABLE_HWLOC:BOOL=ON \
  -D HWLOC_INCLUDE_DIRS:FILEPATH=/home/sems/common/hwloc/current/include \
  -D HWLOC_LIBRARY_DIRS:FILEPATH=/home/sems/common/hwloc/current/lib \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D MPI_BASE_DIR:FILEPATH=/home/sems/common/openmpi/current/ \
  -D CUDA_NVCC_FLAGS:STRING="-arch=sm_20" \
  $EXTRA_ARGS \
  ../Trilinos


