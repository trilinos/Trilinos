/*
* Copyright (c) 2013, Sandia Corporation.
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
* the U.S. Government retains certain rights in this software.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
* 
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
* 
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following
*       disclaimer in the documentation and/or other materials provided
*       with the distribution.
* 
*     * Neither the name of Sandia Corporation nor the names of its
*       contributors may be used to endorse or promote products derived
*       from this software without specific prior written permission.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/ 

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

__global__
void cuda_kernel_saxpy(int num_times, int n, const float* x, float alpha, float* y)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n) {
    for(int j=0; j<num_times; ++j) {
      y[i] = y[i] + alpha*x[i];
    }
  }
}

int cuda_saxpy(int num_times, int n, const float* x, float alpha, float* y)
{
  float *device_x = NULL, *device_y = NULL;
  cudaError_t err = cudaSuccess;

  size_t num_bytes = n*sizeof(float);

  err = cudaMalloc((void**)&device_x, num_bytes);
  if (err != cudaSuccess) {
    std::cerr<<"cudaMalloc(device_x) ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return -1;
  }

  err = cudaMalloc((void**)&device_y, num_bytes);
  if (err != cudaSuccess) {
    std::cerr<<"cudaMalloc(device_y) ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return -1;
  }

  err = cudaMemcpy(device_x, x, num_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    std::cerr<<"cudaMemcpy(device_x, x) ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return -1;
  }

  err = cudaMemcpy(device_y, y, num_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    std::cerr<<"cudaMemcpy(device_y, y) ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return -1;
  }

  //the following two numbers are magic that I got from a cuda-sample program.
  //I need to understand and possibly parameterize these numbers.

  int threadsPerBlock = 256;
  int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
  std::cout<<"Cuda: blocksPerGrid: "<<blocksPerGrid<<", threadsPerBlock: "<<threadsPerBlock<<std::endl;

  cuda_kernel_saxpy<<<blocksPerGrid, threadsPerBlock>>>(num_times, n, device_x, alpha, device_y);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::cerr<<"cuda_kernel_saxpy ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return -1;
  }

  err = cudaMemcpy(y, device_y, num_bytes, cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    std::cerr<<"cudaMemcpy(y, device_y) ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return -1;
  }

  return 0;
}

