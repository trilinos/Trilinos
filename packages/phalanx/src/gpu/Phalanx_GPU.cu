// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Phalanx_GPU.hpp"
#include <iostream>

#include <cuda_runtime_api.h>

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/sort.h"

namespace PHX{

void ListDevices(std::ostream& os)	
{
  int deviceCount = 0;
  cudaGetDeviceCount(&deviceCount);
  for (int dev = 0; dev < deviceCount; ++dev) {
    
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    os << "\nDevice " << dev << ": " << deviceProp.name << std::endl;
	
    os << "  Total amount of global memory: " << deviceProp.totalGlobalMem << std::endl;

    #if CUDART_VERSION >= 2000
    os << "  Number of multiprocessors: " << deviceProp.multiProcessorCount << std::endl;
    //os << "  Number of cores: " << nGpuArchCoresPerSM[deviceProp.major] * deviceProp.multiProcessorCount) << std::endl;
    #endif
    os << "  Total amount of constant memory: " << deviceProp.totalConstMem << std::endl;
    os << "  Total amount of shared memory per block: " << deviceProp.sharedMemPerBlock << std::endl;
    os << "  Total number of registers available per block: " << deviceProp.regsPerBlock << std::endl;
    os << "  Warp size: " << deviceProp.warpSize << std::endl;
    os << "  Maximum number of threads per block: " << deviceProp.maxThreadsPerBlock << std::endl;
    os << "  Maximum sizes of each dimension of a block: " << deviceProp.maxThreadsDim[0] << " " << deviceProp.maxThreadsDim[1] 
       << " " << deviceProp.maxThreadsDim[2] << std::endl;
    os << "  Maximum sizes of each dimension of a grid: " << deviceProp.maxGridSize[0] << " " << deviceProp.maxGridSize[1] 
       << " " << deviceProp.maxGridSize[2] << std::endl;
    os << "" << dev << std::endl;
    os << "" << dev << std::endl;

/*
        shrLog("  Maximum memory pitch:                          %u bytes\n", deviceProp.memPitch);
        shrLog("  Texture alignment:                             %u bytes\n", deviceProp.textureAlignment);
        shrLog("  Clock rate:                                    %.2f GHz\n", deviceProp.clockRate * 1e-6f);
    #if CUDART_VERSION >= 2000
        shrLog("  Concurrent copy and execution:                 %s\n", deviceProp.deviceOverlap ? "Yes" : "No");
    #endif
    #if CUDART_VERSION >= 2020
        shrLog("  Run time limit on kernels:                     %s\n", deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
        shrLog("  Integrated:                                    %s\n", deviceProp.integrated ? "Yes" : "No");
        shrLog("  Support host page-locked memory mapping:       %s\n", deviceProp.canMapHostMemory ? "Yes" : "No");
        shrLog("  Compute mode:                                  %s\n", deviceProp.computeMode == cudaComputeModeDefault ?
			                                                            "Default (multiple host threads can use this device simultaneously)" :
		                                                                deviceProp.computeMode == cudaComputeModeExclusive ?
																		"Exclusive (only one host thread at a time can use this device)" :
		                                                                deviceProp.computeMode == cudaComputeModeProhibited ?
																		"Prohibited (no host thread can use this device)" :
																		"Unknown");
    #endif
    #if CUDART_VERSION >= 3000
        shrLog("  Concurrent kernel execution:                   %s\n", deviceProp.concurrentKernels ? "Yes" : "No");
    #endif
    #if CUDART_VERSION >= 3010
        shrLog("  Device has ECC support enabled:                %s\n", deviceProp.ECCEnabled ? "Yes" : "No");
    #endif
	}

*/

  }
}

void EvaluateGPU()
{
  std::cout << "GPU!" << std::endl;
  int deviceCount = 0;
  cudaGetDeviceCount(&deviceCount);
  std::cout << "device count = " << deviceCount << std::endl;
  cudaDeviceProp deviceProp;
  int dev = 0;	
  cudaGetDeviceProperties(&deviceProp, dev);

  // generate 16M random numbers on the host
  thrust::host_vector<int> h_vec(1 << 24);
  //thrust::generate(h_vec.begin(), h_vec.end(), rand);
  // transfer data to the device
  thrust::device_vector<int> d_vec = h_vec;
  // sort data on the device
  thrust::sort(d_vec.begin(), d_vec.end());
  // transfer data back to host
  thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());


}

}
