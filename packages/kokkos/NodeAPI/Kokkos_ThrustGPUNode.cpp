#include "Kokkos_ThrustGPUNode.hpp"
#include <Teuchos_TestForException.hpp>
#include <iostream>
#include <cuda_runtime.h>

namespace Kokkos {

  ThrustGPUNode::ThrustGPUNode(Teuchos::ParameterList &pl)
  {
    using std::cout;
    using std::cerr;
    using std::endl;

    // get node parameters
    int device = pl.get<int>("Device Number",0);
    int verbose = pl.get<int>("Verbose",0);
    // set device
    int deviceCount; cudaGetDeviceCount(&deviceCount); 
    TEST_FOR_EXCEPTION(deviceCount == 0, std::runtime_error,
        "ThrustGPUNode::ThrustGPUNode(): system has no CUDA devices.");
    if (device < 0 || device >= deviceCount) {
      cerr << "ThrustGPUNode::ThrustGPUNode(): specified device number not valid. Using device 0." << endl;
      device = 0;
    }
    cudaDeviceProp deviceProp; 
    cudaSetDevice(device);
    cudaGetDeviceProperties(&deviceProp, device); 
    // as of CUDA 2.1, device prop contains the following fields
    // char name[256]; 
    // size_t totalGlobalMem, sharedMemPerBlock; 
    // int regsPerBlock, warpSize; 
    // size_t memPitch; 
    // int maxThreadsPerBlock, maxThreadsDim[3], maxGridSize[3]; 
    // size_t totalConstMem; 
    // int major, minor;
    // int clockRate; 
    // size_t textureAlignment; 
    // int deviceOverlap; 
    // int multiProcessorCount; 
    // int kernelExecTimeoutEnabled; 
    if (verbose) {
      cout << "ThrustGPUNode attached to device #" << device << " \"" << deviceProp.name 
        << "\", of compute capability " << deviceProp.major << "." << deviceProp.minor
        << endl;
    }
    totalMem_ = deviceProp.totalGlobalMem;
  } 

  ThrustGPUNode::~ThrustGPUNode() {}

  void ThrustGPUNode::sync() const {
    cudaError err = cudaThreadSynchronize();
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error,
        "Kokkos::ThrustGPUNode::sync(): cudaThreadSynchronize() returned error " << err );
  }

}
