
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
  cudaError_t err = cudaSuccess;

  //the following two numbers are magic that I got from a cuda-sample program.
  //I need to understand and possibly parameterize these numbers.

  int threadsPerBlock = 256;
  int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

  cuda_kernel_saxpy<<<blocksPerGrid, threadsPerBlock>>>(num_times, n, x, alpha, y);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::cerr<<"cuda_kernel_saxpy ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return -1;
  }

  return 0;
}

