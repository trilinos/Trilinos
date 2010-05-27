#include "Kokkos_CUDANodeUtils.hpp"
#include "Kokkos_CUDA_util_inline_runtime.h"
#include <iostream>

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
    cutilSafeCallNoSync( cudaFree(ptr) );
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
    node_->allocSize_ -= allocSize_;
#endif
  }

}
