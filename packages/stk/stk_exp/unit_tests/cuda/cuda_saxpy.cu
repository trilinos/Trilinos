
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

