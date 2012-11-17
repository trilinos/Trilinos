#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

# Directory for KokkosArray

KOKKOSARRAY="../../../array"
KOKKOSEMBED="../.."

source ${KOKKOSARRAY}/src/build_common.sh

# Process command line options and set compilation variables
#
# INC_PATH     : required include paths
# CXX          : C++ compiler with compiler-specific options
# CXX_SOURCES  : C++ files for host
# NVCC         : Cuda compiler with compiler-specific options
# NVCC_SOURCES : Cuda sources

#-----------------------------------------------------------------------------

EXECUTABLE="test.exe"

INC_PATH="${INC_PATH} -I. -I${KOKKOSEMBED}/src"

CXX_SOURCES="${CXX_SOURCES} ./TestHost.cpp ./TestMain.cpp"

#-----------------------------------------------------------------------------

if [ -n "${NVCC}" ] ;
then
  NVCC_SOURCES="${NVCC_SOURCES} ./TestCuda.cpp"

  echo ${NVCC} ${INC_PATH} ${NVCC_SOURCES}

  ${NVCC} ${INC_PATH} ${NVCC_SOURCES}
else
  CXX_SOURCES="${CXX_SOURCES} ./TestCuda.cpp"
fi

#-----------------------------------------------------------------------------

echo ${CXX} ${INC_PATH} -o ${EXECUTABLE} ${CXX_SOURCES} ${LIB}

${CXX} -ftree-vectorize ${INC_PATH} -o ${EXECUTABLE} ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------

