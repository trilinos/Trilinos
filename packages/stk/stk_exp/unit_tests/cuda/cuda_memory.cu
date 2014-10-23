
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

#include <vector>

__global__
void cuda_kernel_init(float* x, int n, float init_value)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i < n) {
    x[i] = init_value;
  }
}

float* cuda_alloc_float(size_t n, const float& init_value)
{
  float *device_mem = NULL;
  cudaError_t err = cudaSuccess;

  size_t num_bytes = n*sizeof(float);

  err = cudaMalloc((void**)&device_mem, num_bytes);
  if (err != cudaSuccess) {
    std::cerr<<"cudaMalloc(device_mem) ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return NULL;
  }

  //the following two numbers are magic that I got from a cuda-sample program.
  //I need to understand and possibly parameterize these numbers.

  int threadsPerBlock = 256;
  int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

  cuda_kernel_init<<<blocksPerGrid, threadsPerBlock>>>(device_mem, n, init_value);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::cerr<<"cuda_kernel_init ERROR: "<<cudaGetErrorString(err)<<std::endl;
    return NULL;
  }

  return device_mem;
}

void copy_cuda_memory_to_host(size_t n, const float* device_mem, std::vector<float>& host_vec)
{
  host_vec.resize(n);

  size_t num_bytes = n*sizeof(float);

  cudaError_t err = cudaMemcpy((void*)&host_vec[0], (void*)device_mem, num_bytes, cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    std::cerr<<"cudaMemcpy(host_vec, device_mem) ERROR: "<<cudaGetErrorString(err)<<std::endl;
  }
}

