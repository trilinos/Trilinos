#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

SRC="../src"

INC_PATH="-I. -I${SRC} -I../TPL"
CXX="g++"
CXXFLAGS="-Wall"

CXX_SOURCES=""
CXX_SOURCES="${CXX_SOURCES} ${SRC}/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC}/Host/KokkosArray_Host_Impl.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC}/Host/KokkosArray_Host_MemorySpace.cpp"
CXX_SOURCES="${CXX_SOURCES} UnitTestMain.cpp TestHost.cpp"
CXX_SOURCES="${CXX_SOURCES} ../TPL/gtest/gtest-all.cc"

#-----------------------------------------------------------------------------

while [ -n "${1}" ] ; do

ARG="${1}"
shift 1

case ${ARG} in
#-------------------------------
#----------- OPTIONS -----------
CUDA | Cuda | cuda ) HAVE_CUDA=1 ;;
HWLOC | hwloc ) HAVE_HWLOC=${1} ; shift 1 ;;
OPT | opt | O3 | -O3 ) OPTFLAGS="-O3" ;;
DBG | dbg | g | -g )   OPTFLAGS="-g -DKOKKOSARRAY_BOUNDS_CHECK" ;;
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

if [ -n "${HAVE_CUDA}" ] ;
then
  NVCC_SOURCES="${SRC}/Cuda/*.cu"
  NVCC_SOURCES="${NVCC_SOURCES} TestCudaFunctions.cu"
  CXX_SOURCES="${CXX_SOURCES} TestCuda.cpp"
  LIB="${LIB} -L/usr/local/cuda/lib64 libCuda.a -lcudart -lcuda -lcusparse"
  nvcc -arch=sm_20 -lib -o libCuda.a ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES}
fi

#-----------------------------------------------------------------------------

if [ -n "${HAVE_HWLOC}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ${SRC}/Host/KokkosArray_Host_hwloc.cpp"
  LIB="${LIB} -L${HAVE_HWLOC}/lib -lhwloc"
  INC_PATH="${INC_PATH} -I${HAVE_HWLOC}/include"
else
  CXX_SOURCES="${CXX_SOURCES} ${SRC}/Host/KokkosArray_Host_hwloc_unavailable.cpp"
fi
#-----------------------------------------------------------------------------
# Option for PTHREAD or WINTHREAD eventually

HAVE_PTHREAD=1

if [ -n "${HAVE_PTHREAD}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ${SRC}/Host/KokkosArray_Host_pthread.cpp"
  LIB="${LIB} -lpthread"
else
  CXX_SOURCES="${CXX_SOURCES} ${SRC}/Host/KokkosArray_Host_nothread.cpp"
fi

#-----------------------------------------------------------------------------

echo "Building regular files as: " ${CXX} ${CXXFLAGS} ${OPTFLAGS}

${CXX} ${CXXFLAGS} ${OPTFLAGS} ${INC_PATH} ${TEST_MACRO} -o unit_test.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a ThreadPool_config.h

#-----------------------------------------------------------------------------


