#include "Kokkos_CUDANodeUtils.hpp"
#include "Kokkos_CUDA_util_inline_runtime.h"

namespace Kokkos {

  void CUDANodeDeallocator::free(void *ptr) {
    cutilSafeCallNoSync( cudaFree(ptr) );
  }

}
