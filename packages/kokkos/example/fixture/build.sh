#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

# Directory for Kokkos

KOKKOS="../../core"

source ${KOKKOS}/src/build_common.sh

# Process command line options and set compilation variables
#
# INC_PATH     : required include paths
# CXX          : C++ compiler with compiler-specific options
# CXX_SOURCES  : C++ files for host
# NVCC         : Cuda compiler with compiler-specific options
# NVCC_SOURCES : Cuda sources

#-----------------------------------------------------------------------------

INC_PATH="${INC_PATH}"

CXX_SOURCES="${CXX_SOURCES} *.cpp"

#-----------------------------------------------------------------------------
# If Cuda compiler set build Cuda source code

if [ -n "${NVCC}" ] ;
then
  NVCC_SOURCES="${NVCC_SOURCES} TestFixture.cu"
  CXX_SOURCES="${CXX_SOURCES}"

  echo ${NVCC} ${INC_PATH} ${NVCC_SOURCES}

  ${NVCC} ${INC_PATH} ${NVCC_SOURCES}
fi

#-----------------------------------------------------------------------------
# Build C++ source code:

EXEC="example.exe"

echo ${CXX} ${INC_PATH} -o ${EXEC} ${CXX_SOURCES} ${LIB}

rm -f testit.exe

${CXX} ${INC_PATH} -o ${EXEC} ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------


