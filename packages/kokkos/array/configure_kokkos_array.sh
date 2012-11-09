#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
# Additional command-line arguments given to this script will be
# passed directly to CMake.
#
CMAKE_EXTRA_ARGS=$@

#
# Force CMake to re-evaluate build options.
#
rm -rf CMakeCache.txt

#-----------------------------------------------------------------------------
# Location of Trilinos source tree.

CMAKE_TRILINOS_BASE_DIR="../Trilinos"

#-----------------------------------------------------------------------------
# MPI configuation:
#
# Must have the MPI_BASE_DIR so that the
# include path can be passed to the Cuda compiler

CMAKE_MPI=""
CMAKE_MPI="${CMAKE_MPI} -D TPL_ENABLE_MPI:BOOL=ON"
CMAKE_MPI="${CMAKE_MPI} -D MPI_BASE_DIR:PATH=/home/sems/common/openmpi/current"

#-----------------------------------------------------------------------------
# Hardware locality cmake configuration:

HWLOC_BASE_DIR="/home/sems/common/hwloc/current"

CMAKE_HWLOC=""
CMAKE_HWLOC="${CMAKE_HWLOC} -D TPL_ENABLE_HWLOC:BOOL=ON"
CMAKE_HWLOC="${CMAKE_HWLOC} -D HWLOC_INCLUDE_DIRS:FILEPATH=${HWLOC_BASE_DIR}/include"
CMAKE_HWLOC="${CMAKE_HWLOC} -D HWLOC_LIBRARY_DIRS:FILEPATH=${HWLOC_BASE_DIR}/lib"

#-----------------------------------------------------------------------------
# Cuda cmake configuration:
#
# Note:  Must turn off CUDA_PROPAGATE_HOST_FLAGS because the
#        Tribits wrapper on cmake forces -pedantic, which results in
#        a flood of warnings from nvcc compiler produced code.
#        This means compiler options must be passed manually.
#
# Note:  Options to CUDA_NVCC_FLAGS must be semi-colon delimited,
#        this is different than the standard CMAKE_CXX_FLAGS syntax.

# Cuda compilation flags:

CUDA_NVCC_FLAGS="-arch=sm_20"
CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-DCUDA_RELEASE_VERSION=`nvcc --version | sed -n -e '/release/{s/^.*release //;s/,.*$//;s/\./0/;p}'`"
CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-Xcompiler;-Wall,-ansi"
CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-O3"

CMAKE_CUDA=""
CMAKE_CUDA="${CMAKE_CUDA} -D TPL_ENABLE_CUDA:BOOL=ON"
CMAKE_CUDA="${CMAKE_CUDA} -D TPL_ENABLE_CUSPARSE:BOOL=ON"
CMAKE_CUDA="${CMAKE_CUDA} -D CUDA_VERBOSE_BUILD:BOOL=OFF"
CMAKE_CUDA="${CMAKE_CUDA} -D CUDA_PROPAGATE_HOST_FLAGS:BOOL=OFF"
CMAKE_CUDA="${CMAKE_CUDA} -D CUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS}"

#-----------------------------------------------------------------------------
# KokkosArray cmake configuration to use Pthreads, HWLOC, and Cuda

CMAKE_KOKKOSARRAY=""
CMAKE_KOKKOSARRAY="${CMAKE_KOKKOSARRAY} ${CMAKE_CUDA}"
CMAKE_KOKKOSARRAY="${CMAKE_KOKKOSARRAY} ${CMAKE_HWLOC}"
CMAKE_KOKKOSARRAY="${CMAKE_KOKKOSARRAY} -D Trilinos_ENABLE_KokkosArray:BOOL=ON"
CMAKE_KOKKOSARRAY="${CMAKE_KOKKOSARRAY} -D KokkosArray_ENABLE_DeviceCuda:BOOL=ON"
CMAKE_KOKKOSARRAY="${CMAKE_KOKKOSARRAY} -D TPL_ENABLE_Pthread:BOOL=ON"

#-----------------------------------------------------------------------------

cmake \
  ${CMAKE_MPI} \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="-Werror" \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
  ${CMAKE_KOKKOSARRAY} \
  ${CMAKE_EXTRA_ARGS} \
  ${CMAKE_TRILINOS_BASE_DIR}

#-----------------------------------------------------------------------------

