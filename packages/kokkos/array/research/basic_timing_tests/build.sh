#!/bin/bash
export INC_PATH="-I. -I../../src"
export CXX_SOURCES="../../src/impl/*.cpp"

case ${1} in
HOST | Host | host )
  HOST_SOURCES="../../src/DeviceHost/*.cpp"
  g++ -Wall -o host_test.exe ${INC_PATH} ${CXX_SOURCES} ${HOST_SOURCES} mainHost.cpp
  ;;
CUDA | Cuda | cuda )
  export NVCC_SOURCES="../../src/DeviceCuda/*.cu"
  nvcc -arch=sm_20 -lcudart -lcusparse -o cuda_test.exe ${INC_PATH} ${NVCC_SOURCES} ${CXX_SOURCES} mainCuda.cu
  ;;
*)
  echo 'unknown options: ' ${1}
esac

