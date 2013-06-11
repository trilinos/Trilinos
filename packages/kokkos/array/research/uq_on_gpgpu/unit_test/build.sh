#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

if [ "STOKHOS" = "${1}" ] ;
then
  shift 1
  HAVE_STOKHOS=${1}
  shift 1
fi

#-----------------------------------------------------------------------------
# Directory for KokkosArray

KOKKOSARRAY="../../.."

source ${KOKKOSARRAY}/src/build_common.sh

# Process command line options and set compilation variables
#
# INC_PATH     : required include paths
# CXX          : C++ compiler with compiler-specific options
# CXX_SOURCES  : C++ files for host
# NVCC         : Cuda compiler with compiler-specific options
# NVCC_SOURCES : Cuda sources

#-----------------------------------------------------------------------------
# Paths and sources

WORK_PATH="../src"
INC_PATH="${INC_PATH} -I. -I${WORK_PATH}"

CXX_SOURCES="${CXX_SOURCES} ${WORK_PATH}/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} main.cpp TestHost.cpp"
CXX_SOURCES="${CXX_SOURCES} TestLegendre.cpp"

#-----------------------------------------------------------------------------
# Put path to Trilinos top-level build directory after "stokhos" on command
# line to enable the original matrix-free algorithm, which uses stokhos
NVCC_FLAGS=""
if [ -n "${HAVE_STOKHOS}" ]
then
    if [ ! -d ${HAVE_STOKHOS} ] ;
    then
      echo "${HAVE_STOKHOS} does not exist"
      exit 1
    fi

    TRILINOS_BUILD_PATH="${HAVE_STOKHOS}"
    #TRILINOS_BUILD_PATH="/home/etphipp/Trilinos/build/opt_serial_cuda"
    
    TEUCHOS_SRC="../../../../../teuchos"
    TEUCHOS_BIN="${TRILINOS_BUILD_PATH}/packages/teuchos"
    TEUCHOS_CORE="-I${TEUCHOS_SRC}/core/src -I${TEUCHOS_BIN}/core/src"
    TEUCHOS_PARM="-I${TEUCHOS_SRC}/parameterlist/src -I${TEUCHOS_BIN}/parameterlist/src"
    TEUCHOS_COMM="-I${TEUCHOS_SRC}/comm/src -I${TEUCHOS_BIN}/comm/src"
    TEUCHOS_NUM="-I${TEUCHOS_SRC}/numerics/src -I${TEUCHOS_BIN}/numerics/src"
    TEUCHOS_INC="${TEUCHOS_CORE} ${TEUCHOS_PARM} ${TEUCHOS_COMM} ${TEUCHOS_NUM}"
    STOKHOS_INC="-I../../../../../stokhos/src -I${TRILINOS_BUILD_PATH}/packages/stokhos/src"
    TEUCHOS_LIB="-L${TEUCHOS_BIN}/core/src -L${TEUCHOS_BIN}/parameterlist/src -L${TEUCHOS_BIN}/comm/src -L${TEUCHOS_BIN}/numerics/src -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore"
    INC_PATH="${INC_PATH} ${TEUCHOS_INC} ${STOKHOS_INC}"
    LAPACK="/usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3 -lm"
    LIB="${LIB} ${TEUCHOS_LIB} ${LAPACK}"
    CXX="${CXX} -DHAVE_KOKKOSARRAY_STOKHOS"
    NVCC_FLAGS="-DHAVE_KOKKOSARRAY_STOKHOS"
fi

#-----------------------------------------------------------------------------
# Option for CUDA

if [ -n "${NVCC}" ] ;
then
  NVCC="${NVCC} --generate-line-info -Xptxas -v"
  NVCC_SOURCES="${NVCC_SOURCES} TestCuda.cu"

  ${NVCC} ${NVCC_FLAGS} ${INC_PATH} ${NVCC_SOURCES}

else
  CXX_SOURCES="${CXX_SOURCES} TestCudaStub.cpp"
fi

#-----------------------------------------------------------------------------

echo ${CXX}

${CXX} ${INC_PATH} -o test_uq.exe ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------



