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

