// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stk_util/util/AlignedAllocator.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <iostream>
#include "Kokkos_Core.hpp"

namespace non_std {
namespace impl {

#ifdef KOKKOS_ENABLE_CUDA
void* CUDAPinnedAndMappedAlignedAllocate(size_t size)
{
  void* ret;
  cudaError_t status = cudaHostAlloc(&ret, size, cudaHostAllocMapped);

  STK_ThrowRequireMsg(status == cudaSuccess, "Error during CUDAPinnedAndMappedAlignedAllocate: " + std::string(cudaGetErrorString(status)));
  return ret;
}

void CUDAPinnedAndMappedAlignedDeallocate(void* p)
{
  cudaError_t status = cudaFreeHost(p);
  STK_ThrowRequireMsg(status == cudaSuccess, "Error during CUDAPinnedAndMappedAlignedDellocate: " + std::string(cudaGetErrorString(status)));
}
#endif

#ifdef KOKKOS_ENABLE_HIP
void* HIPPinnedAndMappedAlignedAllocate(size_t size)
{
  void* ret;
  hipError_t status = hipHostMalloc(&ret, size, hipHostMallocMapped);

  STK_ThrowRequireMsg(status == hipSuccess, "Error during hipPinnedAndMappedAlignedAllocate: " + std::string(hipGetErrorString(status)));
  return ret;
}

void HIPPinnedAndMappedAlignedDeallocate(void* p)
{
  if (!Kokkos::is_initialized())
    std::cerr << "Error during HIPPinnedAndMappedAlignedDeallocate::hipFreeHost: Kokkos not initialized"<< std::endl;

  hipError_t status = hipHostFree(p);
  STK_ThrowRequireMsg(status == hipSuccess, "Error during HIPPinnedAndMappedAlignedDellocate: " + std::string(hipGetErrorString(status)));
}
#endif
}
}
