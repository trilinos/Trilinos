#include "Kokkos_CUDANodeUtils.hpp"
#include <iostream>
#include <cuda_runtime.h>

namespace Kokkos {

  CUDANodeDeallocator::CUDANodeDeallocator(size_t sizeInBytes, const Teuchos::RCP<CUDANodeMemoryModel> &node)
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
  : node_(node)
  , allocSize_(sizeInBytes)
#endif
  {
    (void)sizeInBytes;
    (void)node;
  }

  void CUDANodeDeallocator::free(void *ptr) {
    cudaError_t err = cudaFree(ptr);
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error, 
        "Kokkos::CUDANodeDeallocator::free(): cudaFree() returned error:\n"
        << cudaGetErrorString(err) 
      );
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
    node_->allocSize_ -= allocSize_;
#endif
  }

}
