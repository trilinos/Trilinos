#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

KOKKOSARRAY="../../array"

source ${KOKKOSARRAY}/src/build_common.sh

# Process command line options and set compilation variables
#
# INC_PATH     : required include paths
# CXX          : C++ compiler with compiler-specific options
# CXX_SOURCES  : C++ files for host
# NVCC         : Cuda compiler with compiler-specific options
# NVCC_SOURCES : Cuda sources

#-----------------------------------------------------------------------------

INC_PATH="${INC_PATH} -I. -I../src"

CXX_SOURCES="${CXX_SOURCES} TestArrayExp.cpp TestMain.cpp"

#-----------------------------------------------------------------------------

if [ -n "${NVCC}" ] ;
then
  NVCC_SOURCES="${NVCC_SOURCES} TestArrayExp.cpp"
  NVCC="${NVCC} -DTEST_KOKKOSARRAY_DEVICE=KokkosArray::Cuda"

  echo ${NVCC} ${INC_PATH} ${NVCC_SOURCES}

  ${NVCC} ${INC_PATH} ${NVCC_SOURCES}
fi

#-----------------------------------------------------------------------------

CXX="${CXX} -DTEST_KOKKOSARRAY_DEVICE=KokkosArray::Host"

echo ${CXX} ${INC_PATH} -o unit_test.exe ${CXX_SOURCES} ${LIB}

${CXX} ${INC_PATH} -o unit_test.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------


