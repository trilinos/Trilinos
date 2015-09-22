
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

int check_cuda()
{
  int deviceCount = 0;
  cudaError_t err = cudaGetDeviceCount(&deviceCount);
  std::cout<<"Detected "<<deviceCount<<" cuda devices."<<std::endl;

  if (deviceCount > 0) {
    int deviceId = 0;
    cudaSetDevice(deviceId);
//    cudaDeviceReset();

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, deviceId);

    size_t cudaFreeMem = 0, cudaTotalMem = 0;
    err = cudaMemGetInfo(&cudaFreeMem, &cudaTotalMem);
    if (err != cudaSuccess) {
      std::cerr<<"check_cuda: ERROR, cudaMemGetInfo: "<<cudaGetErrorString(err)<<std::endl;
      return -1;
    }

    double MB = 1024*1024;
    double freeMB = cudaFreeMem/MB;
    double totalMB = cudaTotalMem/MB;

    std::cout<<"cuda device '"<<deviceProp.name<<"', selected properties:\n";
    std::cout<<"\t"<<freeMB<<"MB free out of "<<totalMB<<"MB total memory\n";
    std::cout<<"\t"<<deviceProp.multiProcessorCount<<" multi-processors, maxThreadsPerBlock: "<<deviceProp.maxThreadsPerBlock<<", warpSize: "<<deviceProp.warpSize<<std::endl;
  }

  return deviceCount > 0 ? 0 : -1;//return 0 if device(s) found, -1 if not
}

