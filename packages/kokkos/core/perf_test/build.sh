#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

# Directory for Kokkos

KOKKOS=".."

source ${KOKKOS}/src/build_common.sh

# Process command line options and set compilation variables
#
# INC_PATH     : required include paths
# CXX          : C++ compiler with compiler-specific options
# CXX_SOURCES  : C++ files for host
# NVCC         : Cuda compiler with compiler-specific options
# NVCC_SOURCES : Cuda sources

#-----------------------------------------------------------------------------
# Add TPL and local source code:

CXX_SOURCES="${CXX_SOURCES} ${KOKKOS}/../TPL/gtest/gtest-all.cc"
CXX_SOURCES="${CXX_SOURCES} PerfTestMain.cpp PerfTestHost.cpp"

#-----------------------------------------------------------------------------

if [ -n "${NVCC}" ] ;
then
  NVCC_SOURCES="${NVCC_SOURCES} PerfTestCuda.cu"
  CXX_SOURCES="${CXX_SOURCES} PerfTestCuda.cpp"

  echo ${NVCC} ${INC_PATH} ${NVCC_SOURCES}

  ${NVCC} ${INC_PATH} ${NVCC_SOURCES}
fi

#-----------------------------------------------------------------------------

rm -f perf_test.exe

echo ${CXX} ${INC_PATH} -o perf_test.exe ${CXX_SOURCES} ${LIB}

${CXX} ${INC_PATH} -o perf_test.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------


