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
#-------------------------------
#----------- OPTIONS -----------
OPT | opt | O3 | -O3 )
  CXXFLAGS="${CXXFLAGS} -O3"
  ;;
DBG | dbg | g | -g )
  CXXFLAGS="${CXXFLAGS} -g"
  ;;
#-------------------------------
#---------- COMPILERS ----------
GNU | gnu | g++ )
  export CXX="g++"
  export CXXFLAGS="-Wall"
  ;;
INTEL | intel | icc )
  export CXX="icc"
  export CXXFLAGS="-Wall"
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
  ${NVCC} ${CXXFLAGS} ${INC_PATH} ${NVCC_SOURCES} ;
fi

echo ${CXX} ${CXXFLAGS}

${CXX} ${CXXFLAGS} ${INC_PATH} ${TEST_MACRO} -o mini_test.exe main.cpp ${CXX_SOURCES} ${LIB}

rm -f *.o *.a ThreadPool_config.h

#-----------------------------------------------------------------------------

