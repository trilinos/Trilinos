#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
#-----------------------------------------------------------------------------
# Location of Trilinos source tree.

TRILINOS_SOURCE_DIRECTORY="../Trilinos"

# If installing the location to install
# TRILINOS_INSTALL_DIRECTORY="../TrilinosInstall"

# Build options:

# OPTIMIZATION="DEBUG"
OPTIMIZATION="RELEASE"

# HWLOC_DIRECTORY="/home/sems/common/hwloc/current"

# MPI_DIRECTORY="/home/sems/common/openmpi/current"

# PTHREAD="ON"
# OPENMP="ON"

# INTEL="ON"
# INTEL_XEON_PHI="ON"

# CUDA_ARCH="20"
# CUDA_ARCH="30"
# CUDA_ARCH="35"

# ECLIPSE="ON"

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

CMAKE_CONFIGURE=""
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"

if [ -n "${TRILINOS_INSTALL_DIRECTORY}" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIRECTORY}"
fi

KOKKOS_LINALG_CXXFLAG="-DKOKKOS_FAST_COMPILE"

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_FLAGS:STRING=${KOKKOS_LINALG_FASTCOMPILE}"

#-----------------------------------------------------------------------------
# MPI configuation:
#
# Must have the MPI_BASE_DIR so that the
# include path can be passed to the Cuda compiler

if [ -n "${MPI_DIRECTORY}" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_MPI:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D MPI_BASE_DIR:PATH=${MPI_DIRECTORY}"
fi

#-----------------------------------------------------------------------------
# Pthread configuation:

if [ "${PTHREAD}" = "ON" ] ;
then 
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_Pthread:BOOL=ON"
fi

#-----------------------------------------------------------------------------
# OpenMP configuation:

if [ "${OPENMP}" = "ON" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_OpenMP:BOOL=ON"
fi

#-----------------------------------------------------------------------------
# Hardware locality cmake configuration:

if [ -n "${HWLOC_DIRECTORY}" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_HWLOC:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HWLOC_INCLUDE_DIRS:FILEPATH=${HWLOC_DIRECTORY}/include"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HWLOC_LIBRARY_DIRS:FILEPATH=${HWLOC_DIRECTORY}/lib"
fi

#-----------------------------------------------------------------------------
# Intel compiler:

if [ "${INTEL}" = "ON" -o "${INTEL_XEON_PHI}" = "ON" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_C_COMPILER=icc"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_COMPILER=icpc"
fi

#-----------------------------------------------------------------------------
# Cross-compile for Intel Xeon Phi:

if [ "${INTEL_XEON_PHI}" = "ON" ] ;
then

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_SYSTEM_NAME=Linux"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_FLAGS:STRING=-mmic"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_C_FLAGS:STRING=-mmic"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_Fortran_COMPILER:FILEPATH=ifort"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D BLAS_LIBRARY_DIRS:FILEPATH=${MKLROOT}/lib/mic"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D BLAS_LIBRARY_NAMES='mkl_intel_lp64;mkl_sequential;mkl_core;pthread;m'"

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_CHECKED_STL:BOOL=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=''"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D BUILD_SHARED_LIBS:BOOL=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D DART_TESTING_TIMEOUT:STRING=600"

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D LAPACK_LIBRARY_NAMES=''"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_LAPACK_LIBRARIES=''"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_BinUtils=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_Pthread_LIBRARIES=pthread"

  # Cannot cross-compile fortran compatibility checks on the MIC:
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_Fortran:BOOL=OFF"

  # Tell cmake the answers to compile-and-execute tests
  # to prevent cmake from executing a cross-compiled program.
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HAVE_GCC_ABI_DEMANGLE_EXITCODE=0"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HAVE_TEUCHOS_BLASFLOAT_EXITCODE=0"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D LAPACK_SLAPY2_WORKS_EXITCODE=0"
fi

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

if [ -n "${CUDA_ARCH}" ] ;
then
  CUDA_NVCC_FLAGS="-gencode;arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}"
  CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-Xcompiler;-Wall,-ansi"
  CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};${KOKKOS_LINALG_CXXFLAG}"

  if [ "${OPENMP}" = "ON" ] ;
  then
    CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-fopenmp"
  fi

  if [ "${OPTIMIZATION}" = "DEBUG" ] ;
  then
    CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-g"
  else
    CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-O3"
  fi

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_CUDA:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_CUSPARSE:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CUDA_VERBOSE_BUILD:BOOL=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS}"
fi

#-----------------------------------------------------------------------------
# configure Trilinos to only build Kokkos

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TESTS:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosCore:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosContainers:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosLinAlg:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosExample:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosTask:BOOL=ON"

if [ "${OPTIMIZATION}" = "DEBUG" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_BOUNDS_CHECK:BOOL=ON"
fi

#-----------------------------------------------------------------------------

if [ "${ECLIPSE}" = "ON" ] ;
then
  CMAKE_CONFIGURE="-G\"Eclipse CDT4 - Unix Makefiles\" ${CMAKE_CONFIGURE}"
fi

#
# Force CMake to re-evaluate build options.
#
rm -rf CMake* Trilinos* packages Dart* Testing cmake_install.cmake MakeFile*

cmake ${CMAKE_CONFIGURE} ${TRILINOS_SOURCE_DIRECTORY}

#-----------------------------------------------------------------------------

