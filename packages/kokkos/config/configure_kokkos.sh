#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
# Additional command-line arguments given to this script will be
# passed directly to CMake.
#

#
# Force CMake to re-evaluate build options.
#
rm -rf CMake* Trilinos* packages Dart* Testing cmake_install.cmake MakeFile*

#-----------------------------------------------------------------------------
# Location of Trilinos source tree.

CMAKE_PROJECT_DIR="../Trilinos"
CMAKE_INSTALL_PREFIX="../TrilinosInstall"

#-----------------------------------------------------------------------------
# MPI configuation:
#
# Must have the MPI_BASE_DIR so that the
# include path can be passed to the Cuda compiler

CMAKE_MPI=""
CMAKE_MPI="${CMAKE_MPI} -D TPL_ENABLE_MPI:BOOL=ON"
CMAKE_MPI="${CMAKE_MPI} -D MPI_BASE_DIR:PATH=/home/sems/common/openmpi/current"

#-----------------------------------------------------------------------------
# Pthread configuation:

CMAKE_PTHREAD=""
CMAKE_PTHREAD="${CMAKE_PTHREAD} -D TPL_ENABLE_Pthread:BOOL=ON"

#-----------------------------------------------------------------------------
# OpenMP configuation:

CMAKE_OPENMP=""
CMAKE_OPENMP="${CMAKE_OPENMP} -D Trilinos_ENABLE_OpenMP:BOOL=ON"

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

CUDA_NVCC_FLAGS="-gencode;arch=compute_20,code=sm_20"
# CUDA_NVCC_FLAGS="-gencode;arch=compute_30,code=sm_30"
# CUDA_NVCC_FLAGS="-gencode;arch=compute_35,code=sm_35"
CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-Xcompiler;-Wall,-ansi"
CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-O3"

CMAKE_CUDA=""
CMAKE_CUDA="${CMAKE_CUDA} -D TPL_ENABLE_CUDA:BOOL=ON"
CMAKE_CUDA="${CMAKE_CUDA} -D TPL_ENABLE_CUSPARSE:BOOL=ON"
CMAKE_CUDA="${CMAKE_CUDA} -D CUDA_VERBOSE_BUILD:BOOL=OFF"
CMAKE_CUDA="${CMAKE_CUDA} -D CUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS}"

#-----------------------------------------------------------------------------
# configure Trilinos to only bulid Kokkos

CMAKE_TRILINOS="${CMAKE_TRILINOS} -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
CMAKE_TRILINOS="${CMAKE_TRILINOS} -D Trilinos_ENABLE_Fortran:BOOL=OFF"
CMAKE_TRILINOS="${CMAKE_TRILINOS} -D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
CMAKE_TRILINOS="${CMAKE_TRILINOS} -D Trilinos_ENABLE_TESTS:BOOL=ON"
CMAKE_TRILINOS="${CMAKE_TRILINOS} -D Trilinos_ENABLE_KokkosCore:BOOL=ON"
CMAKE_TRILINOS="${CMAKE_TRILINOS} -D Trilinos_ENABLE_KokkosContainers:BOOL=ON"
CMAKE_TRILINOS="${CMAKE_TRILINOS} -D Trilinos_ENABLE_KokkosExample:BOOL=ON"


#-----------------------------------------------------------------------------
# Kokkos cmake configuration to use MPI, Pthreads, HWLOC, and Cuda

CMAKE_CONFIGURE=""
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING='RELEASE'"
# CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING='DEBUG'"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} ${CMAKE_MPI}"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} ${CMAKE_CUDA}"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} ${CMAKE_HWLOC}"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} ${CMAKE_PTHREAD}"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} ${CMAKE_OPENMP}"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} ${CMAKE_TRILINOS}"


#cmake -G"Eclipse CDT4 - Unix Makefiles" ${CMAKE_CONFIGURE} ${CMAKE_PROJECT_DIR}
cmake ${CMAKE_CONFIGURE} ${CMAKE_PROJECT_DIR}

#-----------------------------------------------------------------------------

