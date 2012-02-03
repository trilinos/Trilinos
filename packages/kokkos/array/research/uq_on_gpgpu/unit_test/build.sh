#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------
# Process command line options:

while [ -n "${1}" ] ; do

ARG="${1}"
shift 1

case ${ARG} in
#-------------------------------
#----------- OPTIONS -----------
CUDA | Cuda | cuda ) HAVE_CUDA=1 ;;
HWLOC | hwloc ) HAVE_HWLOC=${1} ; shift 1 ;;
OPT | opt | O3 | -O3 ) OPTFLAGS="-O3" ;;
DBG | dbg | g | -g )   OPTFLAGS="-g" ;;
#-------------------------------
#---------- COMPILERS ----------
GNU | gnu | g++ )
  CXX="g++"
  CXXFLAGS="-Wall"
  ;;
INTEL | intel | icc )
  CXX="icc"
  # -xW = use SSE and SSE2 instructions
  CXXFLAGS="-Wall -xW"
  LIB="${LIB} -lstdc++"
  ;;
#-------------------------------
*) echo 'unknown option: ' ${ARG} ; exit -1 ;;
esac
done

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Paths and sources

SRC_PATH="../../../src"
INC_PATH="-I. -I../src -I${SRC_PATH}"

CXX_SOURCES="main.cpp TestHost.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/Kokkos_Host_Impl.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/Kokkos_Host_MemoryManager.cpp"

#-----------------------------------------------------------------------------

if [ -n "${HAVE_HWLOC}" ] ;
then
  if [ ! -d ${HAVE_HWLOC} ] ;
  then
    echo "${HAVE_HWLOC} does not exist"
    exit 1
  fi

  HWLOC_LIB_PATH="${HAVE_HWLOC}/lib"

  echo "LD_LIBRARY_PATH must include ${HWLOC_LIB_PATH}"

  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/Kokkos_Host_hwloc.cpp"
  LIB="${LIB} -L${HWLOC_LIB_PATH} -lhwloc"
  INC_PATH="${INC_PATH} -I${HAVE_HWLOC}/include"
else
  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/Kokkos_Host_hwloc_unavailable.cpp"
fi

#-----------------------------------------------------------------------------
# Option for PTHREAD or WINTHREAD eventually

HAVE_PTHREAD=1

if [ -n "${HAVE_PTHREAD}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/Kokkos_Host_pthread.cpp"
  LIB="${LIB} -lpthread"
else
  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/Kokkos_Host_nothread.cpp"
fi

#-----------------------------------------------------------------------------
# Option for CUDA

if [ -n "${HAVE_CUDA}" ] ;
then
  TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_CUDA"
  NVCC_PATH="/usr/local/cuda"
  NVCC_SOURCES="${SRC_PATH}/Cuda/*.cu TestCuda.cu"
  NVCC="${NVCC_PATH}/bin/nvcc -arch=sm_20 -lib -o libCuda.a"
  LIB="${LIB} -L${NVCC_PATH}/lib64 libCuda.a -lcudart -lcuda -lcusparse"

  ${NVCC} ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES} ;
else
  CXX_SOURCES="${CXX_SOURCES} TestCudaStub.cpp"
fi

#-----------------------------------------------------------------------------

${CXX} ${CXXFLAGS} ${OPTFLAGS} ${INC_PATH} ${TEST_MACRO} -o test_uq.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------



