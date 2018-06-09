// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Sacado_DynamicArrayTraits.hpp"

#if 0 && defined(HAVE_SACADO_KOKKOSCORE) && defined(KOKKOS_ENABLE_OPENMP)
namespace Sacado {
  namespace Impl {
    const Kokkos::MemoryPool<Kokkos::OpenMP>* global_sacado_openmp_memory_pool = 0;
  }
}
#endif

#if defined(HAVE_SACADO_KOKKOSCORE) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(__CUDACC__)
namespace Sacado {
  namespace Impl {
    const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_host = 0;
    const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_device = 0;
#ifdef KOKKOS_CUDA_USE_RELOCATABLE_DEVICE_CODE
    __device__ const Kokkos::MemoryPool<Kokkos::Cuda>* global_sacado_cuda_memory_pool_on_device = 0;
#endif
  }
}
#endif
