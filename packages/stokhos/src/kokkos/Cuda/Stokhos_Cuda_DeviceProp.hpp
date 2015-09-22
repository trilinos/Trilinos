// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_DEVICE_PROP_HPP
#define STOKHOS_CUDA_DEVICE_PROP_HPP

#include "Kokkos_Core.hpp"

#include "Teuchos_TestForException.hpp"

#include "cuda_runtime_api.h"

namespace Stokhos {

  // Class encapsulating various device attributes
  class DeviceProp {
  public:

    typedef Kokkos::Cuda::size_type size_type;

    size_type compute_capability_major;
    size_type compute_capability_minor;

    size_type shared_memory_capacity;
    size_type shared_memory_granularity;
    size_type max_shmem_per_block;
    size_type max_threads_per_block;
    size_type max_threads_per_sm;
    size_type max_blocks_per_sm;
    size_type max_warps_per_sm;
    size_type warp_size;
    size_type warp_granularity;
    size_type max_regs_per_sm;
    size_type max_regs_per_block;
    size_type reg_bank_size;

    bool has_shuffle;
    bool has_ldg;

    DeviceProp(int device_id = -1) :
      compute_capability_major(0),
      compute_capability_minor(0),
      shared_memory_capacity(0),
      shared_memory_granularity(0),
      max_shmem_per_block(0),
      max_threads_per_block(0),
      max_threads_per_sm(0),
      max_blocks_per_sm(0),
      max_warps_per_sm(0),
      warp_size(0),
      warp_granularity(0),
      max_regs_per_sm(0),
      max_regs_per_block(0),
      reg_bank_size(0),
      has_shuffle(false),
      has_ldg(false)
    {
      // If device_id is negative, use currently selected device
      if (device_id < 0)
        cudaGetDevice(&device_id);

      // Get compute capability
      int major, minor;
      cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor,
                             device_id);
      cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor,
                             device_id);
      compute_capability_major = major;
      compute_capability_minor = minor;

      // Require compute capability >= 2
      TEUCHOS_TEST_FOR_EXCEPTION(
        compute_capability_major < 2, std::logic_error,
        "Cuda compute capability >= 2 is required!");

      // These come from the CUDA occupancy calculator
      if (compute_capability_major == 3) {
        if (compute_capability_minor >= 7) {
          shared_memory_capacity = 112 * 1024;
          max_shmem_per_block = 48 * 1024;
          max_regs_per_sm = 128 * 1024;
          max_regs_per_block = 64 * 1024;
        }
        else {
          shared_memory_capacity = 48 * 1024;
          max_shmem_per_block = 48 * 1024;
          max_regs_per_sm = 64 * 1024;
          max_regs_per_block = 64 * 1024;
        }
        shared_memory_granularity = 256;
        max_threads_per_block = 1024;
        max_threads_per_sm = 2048;
        max_blocks_per_sm = 16;
        max_warps_per_sm = 64;
        warp_size = 32;
        warp_granularity = 4;
        reg_bank_size = 256;
        has_shuffle = true;
        has_ldg = true;
      }

      else if (compute_capability_major == 2) {
        shared_memory_capacity = 48 * 1024;
        shared_memory_granularity = 64;
        max_shmem_per_block = 48 * 1024;
        max_threads_per_block = 1024;
        max_threads_per_sm = 1536;
        max_blocks_per_sm = 8;
        max_warps_per_sm = 48;
        warp_size = 32;
        warp_granularity = 2;
        max_regs_per_sm = 32 * 1024;
        max_regs_per_block = 32 * 1024;
        reg_bank_size = 64;
        has_shuffle = false;
        has_ldg = false;
      }
    }

    // Returns number of registers per thread used by the given kernel
    template <typename Kernel>
    size_type
    get_kernel_registers(Kernel kernel) {
#ifdef __CUDACC__
      typedef void (*func_ptr_t)();
      func_ptr_t func_ptr = reinterpret_cast<func_ptr_t>(kernel);
      cudaFuncAttributes attrib;
      cudaFuncGetAttributes(&attrib, func_ptr);
      return attrib.numRegs;
#else
      return 0;
#endif
    }

    // Returns number of resident warps per sm for the given kernel
    template <typename Kernel>
    size_type
    get_resident_warps_per_sm(Kernel kernel) {
      const size_type regs_per_thread = get_kernel_registers(kernel);
      const size_type regs_per_warp =
        (warp_size*regs_per_thread + reg_bank_size-1) & ~(reg_bank_size-1);
      const size_type warps_per_sm =
        (max_regs_per_sm/regs_per_warp) & ~(warp_granularity-1);
      return warps_per_sm;
    }
  };

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_DEVICE_PROP_HPP */
