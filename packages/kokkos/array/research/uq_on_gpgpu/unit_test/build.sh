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
OPT | opt | O3 | -O3 ) OPTFLAGS="-g -O3" ;;
DBG | dbg | g | -g )   OPTFLAGS="-g" ;;
STOKHOS | stokhos ) HAVE_STOKHOS=${1} ; shift 1 ;;
#-------------------------------
#---------- COMPILERS ----------
GNU | gnu | g++ )
  CXX="g++"
  #CXXFLAGS="-Wall"
  CXXFLAGS=""
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

if [ ! -n "${CXX}" ]
then
  CXX="g++"
  CXXFLAGS="-Wall"
fi

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Paths and sources
WORK_PATH="../src"
SRC_PATH="../../../src"
INC_PATH="-I. -I${WORK_PATH} -I${SRC_PATH}"

# Put path to Trilinos top-level build directory after "stokhos" on command
# line to enable the original matrix-free algorithm, which uses stokhos
if [ -n "${HAVE_STOKHOS}" ]
then
    if [ ! -d ${HAVE_STOKHOS} ] ;
    then
      echo "${HAVE_STOKHOS} does not exist"
      exit 1
    fi

    TRILINOS_BUILD_PATH="${HAVE_STOKHOS}"
    #TRILINOS_BUILD_PATH="/home/etphipp/Trilinos/build/opt_serial_cuda"
    
    TEUCHOS_INC="-I../../../../../teuchos/src -I${TRILINOS_BUILD_PATH}/packages/teuchos/src"
    STOKHOS_INC="-I../../../../../stokhos/src -I${TRILINOS_BUILD_PATH}/packages/stokhos/src"
    TEUCHOS_LIB="${TRILINOS_BUILD_PATH}/packages/teuchos/src"
    INC_PATH="${INC_PATH} ${TEUCHOS_INC} ${STOKHOS_INC}"
    LIB="${LIB} -L${TEUCHOS_LIB} -lteuchos /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3 -lm"
    CXXFLAGS="${CXXFLAGS} -DHAVE_KOKKOSARRAY_STOKHOS"
fi

CXX_SOURCES="main.cpp TestHost.cpp"
CXX_SOURCES="${CXX_SOURCES} ${WORK_PATH}/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/KokkosArray_Host_Impl.cpp"
CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/KokkosArray_Host_MemorySpace.cpp"

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

  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/KokkosArray_Host_hwloc.cpp"
  LIB="${LIB} -L${HWLOC_LIB_PATH} -lhwloc"
  INC_PATH="${INC_PATH} -I${HAVE_HWLOC}/include"
else
  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/KokkosArray_Host_hwloc_unavailable.cpp"
fi

#-----------------------------------------------------------------------------
# Option for PTHREAD or WINTHREAD eventually

HAVE_PTHREAD=1

if [ -n "${HAVE_PTHREAD}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/KokkosArray_Host_pthread.cpp"
  LIB="${LIB} -lpthread"
else
  CXX_SOURCES="${CXX_SOURCES} ${SRC_PATH}/Host/KokkosArray_Host_nothread.cpp"
fi

#-----------------------------------------------------------------------------
# Option for CUDA

if [ -n "${HAVE_CUDA}" ] ;
then
  TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOSARRAY_CUDA"
  NVCC_PATH="/usr/local/cuda"
  NVCC_SOURCES="${SRC_PATH}/Cuda/*.cu TestCuda.cu"
  NVCC="${NVCC_PATH}/bin/nvcc -arch=sm_20 -lib -o libCuda.a"
  LIB="${LIB} -L${NVCC_PATH}/lib64 libCuda.a -lcudart -lcuda -lcusparse"

  ${NVCC} ${CXXFLAGS} ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES} ;
else
  CXX_SOURCES="${CXX_SOURCES} TestCudaStub.cpp"
fi

#-----------------------------------------------------------------------------

${CXX} ${CXXFLAGS} ${OPTFLAGS} ${INC_PATH} ${TEST_MACRO} -o test_uq.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------



