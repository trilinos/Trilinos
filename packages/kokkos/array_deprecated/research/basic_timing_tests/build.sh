#!/bin/bash
export INC_PATH="-I. -I../../src"
export CXX_SOURCES="../../src/impl/*.cpp"

case ${1} in
dirac-host )
  module swap gcc intel
  HOST_SOURCES="../../src/DeviceHost/*.cpp"
  icpc -O3 -o host_test.exe ${INC_PATH} ${CXX_SOURCES} ${HOST_SOURCES} mainHost.cpp
  ;;
dirac-device )
  # Need gcc (version < 4.5)
  module unload gcc/4.5.2
  module load gcc/4.4.2
  export NVCC_SOURCES="../../src/DeviceCuda/*.cu"
  nvcc -arch=sm_20 -lcudart -lcusparse -o cuda_test.exe ${INC_PATH} ${NVCC_SOURCES} ${CXX_SOURCES} mainCuda.cu
  ;;
HOST | Host | host )
  HOST_SOURCES="../../src/DeviceHost/*.cpp"
  g++ -O3 -Wall -o host_test.exe ${INC_PATH} ${CXX_SOURCES} ${HOST_SOURCES} mainHost.cpp
  ;;
CUDA | Cuda | cuda )
  export NVCC_SOURCES="../../src/DeviceCuda/*.cu"
  nvcc -arch=sm_20 -lcudart -lcusparse -o cuda_test.exe ${INC_PATH} ${NVCC_SOURCES} ${CXX_SOURCES} mainCuda.cu
  ;;
*)
  echo 'unknown options: ' ${1}
esac

