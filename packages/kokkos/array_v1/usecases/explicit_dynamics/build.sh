#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

export INC_PATH="-I. -I../../src"
export CXX_SOURCES="../../src/impl/*.cpp ../../src/DeviceHost/*.cpp"
export CXX="g++"
export CXXFLAGS="-Wall"

#-----------------------------------------------------------------------------

while [ -n "${1}" ] ; do

export ARG="${1}"
shift 1

case ${ARG} in
#-------------------------------
#----------- DEVICES -----------
HOST | Host | host )
  export TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_HOST"
  export CXX_SOURCES="${CXX_SOURCES} testHost.cpp"
  ;;
PTHREAD | Pthread | pthread )
  export TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_PTHREAD"
  export CXX_SOURCES="${CXX_SOURCES} ../../src/DevicePthread/*.cpp testPthread.cpp"
  export LIB="${LIB} -lpthread"
  ;;
TPI | tpi )
  echo "#define HAVE_PTHREAD" > ThreadPool_config.h
  export TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_TPI"
  export CXX_SOURCES="${CXX_SOURCES} ../../src/DeviceTPI/*.cpp testTPI.cpp"
  export CXX_SOURCES="${CXX_SOURCES} ../../../../ThreadPool/src/*.c"
  export INC_PATH="${INC_PATH} -I../../../../ThreadPool/src"
  export LIB="${LIB} -lpthread"
  ;;
CUDA | Cuda | cuda )
  export TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_CUDA"
  export NVCC_SOURCES="../../src/DeviceCuda/*.cu testCuda.cu"
  export NVCC_PATH="/usr/local/cuda"
  export NVCC="${NVCC_PATH}/bin/nvcc -arch=sm_20 -lib -o libCuda.a"
  export LIB="${LIB} -L${NVCC_PATH}/lib64 libCuda.a -lcudart -lcuda -lcusparse"
  ;;
TBB | tbb )
  export TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_TBB"
  export CXX_SOURCES="${CXX_SOURCES} ../../src/DeviceTBB/*.cpp testTBB.cpp"
  export LIB="${LIB} -ltbb"
  ;;
NUMA | numa )
  export TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_NUMA"
  export CXX_SOURCES="${CXX_SOURCES} testNUMA.cpp"
  export CXX_SOURCES="${CXX_SOURCES} ../../src/DeviceNUMA/Kokkos_DeviceNUMA.cpp"
  export CXX_SOURCES="${CXX_SOURCES} ../../src/DeviceNUMA/Kokkos_DeviceNUMA_hwloc.cpp"
  export CXX_SOURCES="${CXX_SOURCES} ../../src/DeviceNUMA/Kokkos_DeviceNUMA_pthread.cpp"
  export INC_PATH="${INC_PATH} -I${1}/include"
  export LIB="${LIB} -L${1}/lib -lhwloc -lpthread"
  shift 1
  ;;
#-------------------------------
#----------- OPTIONS -----------
OPT | opt | O3 | -O3 )
  CXXFLAGS="${CXXFLAGS} -O3"
  NVCCFLAGS="${NVCCFLAGS} -O3"
  ;;
DBG | dbg | g | -g )
  CXXFLAGS="${CXXFLAGS} -g"
  NVCCFLAGS="${NVCCFLAGS} -g"
  ;;
#-------------------------------
#---------- COMPILERS ----------
GNU | gnu | g++ )
  export CXX="g++"
  export CXXFLAGS="-Wall"
  ;;
INTEL | intel | icc )
  export CXX="icc"
  # -xW = use SSE and SSE2 instructions
  export CXXFLAGS="-Wall -xW"
  export LIB="${LIB} -lstdc++"
  ;;
#-------------------------------
*)
  echo 'unknown option: ' ${ARG}
  exit -1
esac
done

#-----------------------------------------------------------------------------

if [ -n "${NVCC}" ] ;
then
  echo "Building CUDA files as: " ${NVCC} ${NVCCFLAGS}
  ${NVCC} ${NVCCFLAGS} ${INC_PATH} ${NVCC_SOURCES} ;
fi

echo "Building regular files as: " ${CXX} ${CXXFLAGS} ${INC_PATH} ${LIB}

${CXX} ${CXXFLAGS} ${INC_PATH} ${TEST_MACRO} -o explicit_dynamics.x main.cpp ${CXX_SOURCES} ${LIB}

rm -f *.o *.a ThreadPool_config.h
