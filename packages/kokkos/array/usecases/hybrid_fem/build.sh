#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

EXECUTABLE="proxyapp.exe"

INC_PATH="-I. -I../common -I../../src"
CXX="g++"

CXX_SOURCES="./*.cpp"
CXX_SOURCES="${CXX_SOURCES} ../common/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ../../src/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_Impl.cpp"
CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_MemorySpace.cpp"

#-----------------------------------------------------------------------------

while [ -n "${1}" ] ; do

ARG="${1}"
shift 1

case ${ARG} in
#-------------------------------
#----------- OPTIONS -----------
OPT | opt | O3 | -O3 ) OPTFLAGS="${OPTFLAGS} -O3" ;;
DBG | dbg | g | -g )   OPTFLAGS="${OPTFLAGS} -g -DKOKKOSARRAY_BOUNDS_CHECK" ;;
HWLOC | hwloc ) HAVE_HWLOC=${1} ; shift 1 ;;
MPI | mpi )
  HAVE_MPI=${1} ; shift 1
  CXX="${HAVE_MPI}/bin/mpiCC"
  INC_PATH="${INC_PATH} -I${HAVE_MPI}/include"
  OPTFLAGS="${OPTFLAGS} -DHAVE_MPI"
  ;;
CUDA | Cuda | cuda ) HAVE_CUDA=1 ;;
curie )
  CXX="CC"
  HAVE_MPI="/opt/cray/mpt/5.4.4/xt/gemini/mpich2-cray/73"
  INC_PATH="${INC_PATH} -I${HAVE_MPI}/include"
  OPTFLAGS="${OPTFLAGS} -DHAVE_MPI"
  ;;  
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

if [ -n "${HAVE_HWLOC}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_hwloc.cpp"
  LIB="${LIB} -L${HAVE_HWLOC}/lib -lhwloc"
  INC_PATH="${INC_PATH} -I${HAVE_HWLOC}/include"
else
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_hwloc_unavailable.cpp"
fi

#-----------------------------------------------------------------------------
# Option for PTHREAD or WINTHREAD eventually

HAVE_PTHREAD=1

if [ -n "${HAVE_PTHREAD}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_pthread.cpp"
  LIB="${LIB} -lpthread"
else
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_nothread.cpp"
fi

#-----------------------------------------------------------------------------

if [ -n "${HAVE_CUDA}" ] ;
then
  OPTFLAGS="${OPTFLAGS} -DHAVE_CUDA"
  NVCC_SOURCES="./*.cu"
  NVCC_SOURCES="${NVCC_SOURCES} ../../src/Cuda/*.cu"
  LIB="${LIB} -L/usr/local/cuda/lib64 libCuda.a -lcudart -lcuda -lcusparse"
  echo "Building cuda files as: nvcc -arch=sm_20 -lib -o libCuda.a ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES}"
  nvcc -arch=sm_20 -lib -o libCuda.a ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES}
fi

#-----------------------------------------------------------------------------

echo "Building regular files as: " ${CXX} ${CXXFLAGS} ${OPTFLAGS}

${CXX} ${CXXFLAGS} ${OPTFLAGS} ${INC_PATH} -o ${EXECUTABLE} ${CXX_SOURCES} ${LIB}

rm -f *.o *.a

#-----------------------------------------------------------------------------

