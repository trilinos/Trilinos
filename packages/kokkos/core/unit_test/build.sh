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
CXX_SOURCES="${CXX_SOURCES} UnitTestMain.cpp TestThreads.cpp"
CXX_SOURCES="${CXX_SOURCES} TestHWLOC.cpp"

if [ -n "${KOKKOS_HAVE_SERIAL}" ]
then
CXX_SOURCES="${CXX_SOURCES} TestSerial.cpp"
fi

if [ -n "${KOKKOS_HAVE_OPENMP}" ]
then
CXX_SOURCES="${CXX_SOURCES} TestOpenMP.cpp"
fi

#-----------------------------------------------------------------------------
# If Cuda compiler set build Cuda source code

if [ -n "${NVCC}" ] ;
then
  NVCC_SOURCES="${NVCC_SOURCES} TestCudaFunctions.cu"
  CXX_SOURCES="${CXX_SOURCES} TestCuda.cpp"

  echo ${NVCC} ${INC_PATH} ${NVCC_SOURCES}

  ${NVCC} ${INC_PATH} ${NVCC_SOURCES}
fi

#-----------------------------------------------------------------------------
# Build C++ source code:

echo ${CXX} ${INC_PATH} -o unit_test.exe ${CXX_SOURCES} ${LIB}

rm -f unit_test.exe

${CXX} ${INC_PATH} -o unit_test.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------


