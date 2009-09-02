// #include <CudaNodeImpl.hpp>
#include <Kokkos_CUDANode.hpp>
#include <Kokkos_CUDANodeImpl.hpp>
#include <stdexcept>
#include <iostream>

// some CUDA rules of thumb employed here (stolen from slides by Mike Bailey, Oregon State)
// -The number of Blocks should be at least twice the number of MPs 
// -The number of Threads per Block should be a multiple of 64 
// -  192 or 256 are good numbers for Threads/Block 
// We will enforce that numThreads is a power of two (to ease the reduction kernel)
// greater than 64

CUDANode::CUDANode(int device, int numBlocks, int numThreads, int verbose)
: numBlocks_(numBlocks)
, numThreads_(numThreads) 
{
  using std::cout;
  using std::endl;
  using std::runtime_error;
  // enforce that numThreads_ is a multiple of 64
  if (numThreads_ != 64 && numThreads_ != 128 && numThreads_ != 256 && numThreads_ != 512
      && numThreads_ != 1 && numThreads_ != 2 && numThreads_ != 4 && numThreads_ != 8 && numThreads_ != 16
      && numThreads_ != 32) {
    throw runtime_error("CUDANode::CUDANode(): number of threads per block must be a power of two in [1,512].");
  }
  int deviceCount; cudaGetDeviceCount(&deviceCount); 
  if (device >= deviceCount) {
    if (deviceCount == 0) {
      throw runtime_error("CUDANode::CUDANode(): system has no CUDA devices.");
    }
    if (verbose) {
      cout << "CUDANode::CUDANode(): specified device number not valid. Using device 0." << endl;
    }
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
    cout << "CUDANode attached to device #" << device << " \"" << deviceProp.name 
         << "\", of compute capability " << deviceProp.major << "." << deviceProp.minor
         << endl;
  }
  totalMem_ = deviceProp.totalGlobalMem;
} 

CUDANode::~CUDANode() {}

void CUDANode::readyBuffers(CUDABuffer<const void> * const cBuffers,  unsigned int numConstBuffers,
                            CUDABuffer<      void> * const ncBuffers, unsigned int numNonConstBuffers)
{
  for (int i=0; i < numConstBuffers; ++i) {
    if (cBuffers[i].length_ != 0 && cBuffers[i].devc_ptr_ == NULL) {
      throw std::runtime_error("CUDANode::readyBuffers(): not all buffers are ready.");
    }
    else if ( (*cBuffers[i].flags_) == DIRTY_HOST) {
      // char is the same size as void
      CUDABuffer<char> *buff = reinterpret_cast<CUDABuffer<char> *>( cBuffers+i );
      std::cout << "readyBuffers() about to copy buffer to device and mark CLEAN" << std::endl;
      copyToBuffer(buff->length_, buff->host_ptr_, *buff, 0);
      (*buff->flags_) = CLEAN;
    }
  }
  for (int i=0; i < numNonConstBuffers; ++i) {
    if (ncBuffers[i].length_ != 0 && ncBuffers[i].devc_ptr_ == NULL) {
      throw std::runtime_error("CUDANode::readyBuffers(): not all buffers are ready.");
    }
    else if ( (*ncBuffers[i].flags_) == DIRTY_HOST) {
      // char is the same size as void
      CUDABuffer<char> *buff = reinterpret_cast<CUDABuffer<char> *>( ncBuffers+i );
      std::cout << "readyBuffers() about to copy buffer to device and mark DIRTY_DEVC" << std::endl;
      copyToBuffer(buff->length_, buff->host_ptr_, *buff, 0);
      (*buff->flags_) = CLEAN;
    }
    (*ncBuffers[i].flags_) = DIRTY_DEVC;
  }
}

std::map<const void *, CUDABuffer<const void> > CUDANode::out_views_;
