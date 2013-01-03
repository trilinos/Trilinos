#!/bin/bash

#-----------------------------------------------------------------------------
# Shared portion of build script for the base KokkosArray functionality
# Simple build script with options
#-----------------------------------------------------------------------------

if [ ! -d "${KOKKOSARRAY}" ] ;
then
echo "Must set KOKKOSARRAY to the top level KokkosArray directory"
exit -1
fi

#-----------------------------------------------------------------------------

while [ -n "${1}" ] ; do

ARG="${1}"
shift 1

case ${ARG} in
#----------- OPTIONS -----------
OPT | opt | O3 | -O3 ) OPTFLAGS="${OPTFLAGS} -O3" ;;
#-------------------------------
DBG | dbg | g | -g )   OPTFLAGS="${OPTFLAGS} -g -DKOKKOSARRAY_EXPRESSION_CHECK" ;;
#-------------------------------
HWLOC | hwloc ) HAVE_HWLOC=${1} ; shift 1 ;;
#-------------------------------
MPI | mpi )
  HAVE_MPI=${1} ; shift 1
  CXX="${HAVE_MPI}/bin/mpiCC"
  INC_PATH="${INC_PATH} -I${HAVE_MPI}/include"
  OPTFLAGS="${OPTFLAGS} -DHAVE_MPI"
  ;;
#-------------------------------
CUDA | Cuda | cuda )
  HAVE_CUDA="-DHAVE_CUDA"
  OPTFLAGS="${OPTFLAGS} ${HAVE_CUDA}"
  NVCC_SOURCES="${NVCC_SOURCES} ${KOKKOSARRAY}/src/Cuda/*.cu"
  #
  # Extract release major and minor version from compiler
  #
  CUDA_VERSION="`nvcc --version | sed -n -e '/release/{s/^.*release //;s/,.*$//;p}'`"
  CUDA_VERSION_MAJOR=`echo ${CUDA_VERSION} | sed 's/\..*//'`
  CUDA_VERSION_MINOR=`echo ${CUDA_VERSION} | sed 's/^.*\.//'`
  #
  # -x cu : process all files through the Cuda compiler as Cuda code.
  # -lib -o : produce library
  #
  NVCC="nvcc -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30"
  NVCC="${NVCC} -Xcompiler -Wall,-ansi"
  NVCC="${NVCC} -DCUDA_VERSION_MAJOR=${CUDA_VERSION_MAJOR}"
  NVCC="${NVCC} -DCUDA_VERSION_MINOR=${CUDA_VERSION_MINOR}"
  NVCC="${NVCC} -lib -o libCuda.a -x cu"

  LIB="${LIB} libCuda.a -L/usr/local/cuda/lib64 -lcudart -lcuda -lcusparse"
  ;;
#-------------------------------
GNU | gnu | g++ )
  # Turn on lots of warnings and ansi compliance.
  # The Trilinos build system requires '-pedantic'
  # 
  CXX="g++ -Wall -Wextra -ansi -pedantic"
  CXX="${CXX} -rdynamic -DENABLE_TRACEBACK"
  LIB="${LIB} -ldl"
  ;;
#-------------------------------
INTEL | intel | icc )
  # -xW = use SSE and SSE2 instructions
  CXX="icc -Wall -xW"
  LIB="${LIB} -lstdc++"
  ;;
#-------------------------------
MIC | mic )
  CXX="icpc -mmic -ansi-alias -Wall"
  # CXX="${CXX} -mGLOB_default_function_attrs=knc_stream_store_controls=2"
  # CXX="${CXX} -vec-report6"
  # CXX="${CXX} -guide-vec"
  LIB="${LIB} -lstdc++"
  ;;
#-------------------------------
curie )
  CXX="CC"
  HAVE_MPI="/opt/cray/mpt/5.4.4/xt/gemini/mpich2-cray/73"
  INC_PATH="${INC_PATH} -I${HAVE_MPI}/include"
  OPTFLAGS="${OPTFLAGS} -DHAVE_MPI"
  ;;  
#-------------------------------
*) echo 'unknown option: ' ${ARG} ; exit -1 ;;
esac
done

#-----------------------------------------------------------------------------

if [ -z "${CXX}" ] ;
then
  echo "No C++ compiler selected"
  exit -1
fi

CXX="${CXX} ${OPTFLAGS}"

if [ -n "${NVCC}" ] ;
then
  NVCC="${NVCC} ${OPTFLAGS}"
fi

#-----------------------------------------------------------------------------

INC_PATH="${INC_PATH} -I${KOKKOSARRAY}/src"

CXX_SOURCES="${CXX_SOURCES} ${KOKKOSARRAY}/src/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ${KOKKOSARRAY}/src/Host/KokkosArray_Host_Impl.cpp"
CXX_SOURCES="${CXX_SOURCES} ${KOKKOSARRAY}/src/Host/KokkosArray_HostSpace.cpp"

#-----------------------------------------------------------------------------

if [ -n "${HAVE_HWLOC}" ] ;
then

  if [ ! -d ${HAVE_HWLOC} ] ;
  then
    echo "${HAVE_HWLOC} does not exist"
    exit 1
  fi

  echo "LD_LIBRARY_PATH must include ${HAVE_HWLOC}/lib"

  CXX_SOURCES="${CXX_SOURCES} ${KOKKOSARRAY}/src/Host/KokkosArray_Host_hwloc.cpp"
  LIB="${LIB} -L${HAVE_HWLOC}/lib -lhwloc"
  INC_PATH="${INC_PATH} -I${HAVE_HWLOC}/include"
else
  CXX_SOURCES="${CXX_SOURCES} ${KOKKOSARRAY}/src/Host/KokkosArray_Host_hwloc_unavailable.cpp"
fi

#-----------------------------------------------------------------------------
# Option for PTHREAD or WINTHREAD eventually

HAVE_PTHREAD=1

if [ -n "${HAVE_PTHREAD}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ${KOKKOSARRAY}/src/Host/KokkosArray_Host_pthread.cpp"
  LIB="${LIB} -lpthread"
else
  CXX_SOURCES="${CXX_SOURCES} ${KOKKOSARRAY}/src/Host/KokkosArray_Host_nothread.cpp"
fi

#-----------------------------------------------------------------------------

