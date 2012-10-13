#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

# Directory for KokkosArray

KOKKOSARRAY=".."

source ${KOKKOSARRAY}/src/build_common.sh

# Process command line options and set compilation variables
#
# INC_PATH     : required include paths
# CXX          : C++ compiler with compiler-specific options
# CXX_SOURCES  : C++ files for host
# NVCC         : Cuda compiler with compiler-specific options
# NVCC_SOURCES : Cuda sources

#-----------------------------------------------------------------------------
# Add TPL and local source code:

TPL_PATH="${KOKKOSARRAY}/TPL"
INC_PATH="${INC_PATH} -I. -I${TPL_PATH}"

CXX_SOURCES="${CXX_SOURCES} ${TPL_PATH}/gtest/gtest-all.cc"
CXX_SOURCES="${CXX_SOURCES} PerfTestMain.cpp PerfTestHost.cpp"

#-----------------------------------------------------------------------------

if [ -n "${NVCC}" ] ;
then
  NVCC_SOURCES="${NVCC_SOURCES} PerfTestCuda.cu"
  CXX_SOURCES="${CXX_SOURCES} PerfTestCuda.cpp"

  echo ${NVCC} ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES}

  ${NVCC} ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES}
fi

#-----------------------------------------------------------------------------

echo ${CXX} ${OPTFLAGS} ${INC_PATH} -o perf_test.exe ${CXX_SOURCES} ${LIB}

${CXX} ${OPTFLAGS} ${INC_PATH} -o perf_test.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------


