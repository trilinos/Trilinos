#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

# Directory for KokkosArray

KOKKOSARRAY="../.."

source ${KOKKOSARRAY}/src/build_common.sh

# Process command line options and set compilation variables
#
# INC_PATH     : required include paths
# CXX          : C++ compiler with compiler-specific options
# CXX_SOURCES  : C++ files for host
# NVCC         : Cuda compiler with compiler-specific options
# NVCC_SOURCES : Cuda sources

CXX_OPTIONS=""

#-----------------------------------------------------------------------------
# Add local source code:

EXECUTABLE="test_gramschmidt.exe"

INC_PATH="${INC_PATH} -I. -I../common"

CXX_SOURCES="${CXX_SOURCES} ./GramSchmidt.cpp ./Main.cpp"

#-----------------------------------------------------------------------------

if [ -n "${NVCC}" ] ;
then
  NVCC_SOURCES="${NVCC_SOURCES} ./GramSchmidt.cpp"

  echo ${NVCC} ${INC_PATH} ${NVCC_SOURCES}

  ${NVCC} ${INC_PATH} ${NVCC_SOURCES}

  CXX_OPTIONS="-DHAVE_CUDA"
fi

#-----------------------------------------------------------------------------

echo ${CXX} ${CXX_OPTIONS} ${INC_PATH} -o ${EXECUTABLE} ${CXX_SOURCES} ${LIB}

${CXX} ${CXX_OPTIONS} ${INC_PATH} -o ${EXECUTABLE} ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------

