// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_DynamicArrayTraits.hpp"

#if 0 && defined(HAVE_SACADO_KOKKOS) && defined(KOKKOS_ENABLE_OPENMP)
namespace Sacado {
  namespace Impl {
    const Kokkos::MemoryPool<Kokkos::OpenMP>* global_sacado_openmp_memory_pool = 0;
  }
}
#endif

#if defined(HAVE_SACADO_KOKKOS) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(KOKKOS_ENABLE_CUDA) && defined(__CUDACC__)
namespace Sacado {
  namespace Impl {
    const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_host = 0;
    const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_device = 0;
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
    __device__ const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_on_device = 0;
#endif
  }
}
#endif
