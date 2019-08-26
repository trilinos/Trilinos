// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "Tpetra_Details_StaticView.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

namespace { // (anonymous)

// Motivating use cases for the initial size:
//
// 1. GMRES (need restart length (default 30) number of rows)
// 2. Single reduce CG (need 2 x 2)
constexpr size_t minimum_initial_size = sizeof (double) * 30 * 2;

// Intel 17 seems a bit buggy with respect to initialization of
// templated static classes, so let's make the compiler's job really
// easy by having nontemplated static raw pointers.

#ifdef KOKKOS_ENABLE_CUDA

void* cuda_memory_ = nullptr;
size_t cuda_memory_size_ = 0;

void finalize_cuda_memory ()
{
  if (cuda_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::CudaSpace> (cuda_memory_);
    cuda_memory_ = nullptr;
    cuda_memory_size_ = 0;
  }
}

void* cuda_uvm_memory_ = nullptr;
size_t cuda_uvm_memory_size_ = 0;

void finalize_cuda_uvm_memory ()
{
  if (cuda_uvm_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::CudaUVMSpace> (cuda_uvm_memory_);
    cuda_uvm_memory_ = nullptr;
    cuda_uvm_memory_size_ = 0;
  }
}

void* cuda_host_pinned_memory_ = nullptr;
size_t cuda_host_pinned_memory_size_ = 0;

void finalize_cuda_host_pinned_memory ()
{
  if (cuda_host_pinned_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::CudaHostPinnedSpace> (cuda_host_pinned_memory_);
    cuda_host_pinned_memory_ = nullptr;
    cuda_host_pinned_memory_size_ = 0;
  }
}
#endif // KOKKOS_ENABLE_CUDA

void* host_memory_ = nullptr;
size_t host_memory_size_ = 0;

void finalize_host_memory ()
{
  if (host_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::HostSpace> (host_memory_);
    host_memory_ = nullptr;
    host_memory_size_ = 0;
  }
}

} // namespace (anonymous)

#ifdef KOKKOS_ENABLE_CUDA

void*
StaticKokkosAllocation<Kokkos::CudaSpace>::
resize (Kokkos::CudaSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::CudaSpace;
  static bool created_finalize_hook = false;

  if (size > cuda_memory_size_) {
    if (cuda_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (cuda_memory_);
    }
    const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
    cuda_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    cuda_memory_size_ = size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_cuda_memory);
    created_finalize_hook = true;
  }

  return cuda_memory_;
}

void*
StaticKokkosAllocation<Kokkos::CudaUVMSpace>::
resize (Kokkos::CudaUVMSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::CudaUVMSpace;
  static bool created_finalize_hook = false;

  const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
  if (req_size > cuda_uvm_memory_size_) {
    if (cuda_uvm_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (cuda_uvm_memory_);
    }
    cuda_uvm_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    cuda_uvm_memory_size_ = req_size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_cuda_uvm_memory);
    created_finalize_hook = true;
  }

  return cuda_uvm_memory_;
}

void*
StaticKokkosAllocation<Kokkos::CudaHostPinnedSpace>::
resize (Kokkos::CudaHostPinnedSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::CudaHostPinnedSpace;
  static bool created_finalize_hook = false;

  const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
  if (req_size > cuda_host_pinned_memory_size_) {
    if (cuda_host_pinned_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (cuda_host_pinned_memory_);
    }
    cuda_host_pinned_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    cuda_host_pinned_memory_size_ = req_size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_cuda_host_pinned_memory);
    created_finalize_hook = true;
  }

  return cuda_host_pinned_memory_;
}

#endif // KOKKOS_ENABLE_CUDA

void*
StaticKokkosAllocation<Kokkos::HostSpace>::
resize (Kokkos::HostSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::HostSpace;
  static bool created_finalize_hook = false;

  const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
  if (req_size > host_memory_size_) {
    if (host_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (host_memory_);
    }
    host_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    host_memory_size_ = req_size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_host_memory);
    created_finalize_hook = true;
  }

  return host_memory_;
}

} // namespace Impl
} // namespace Details
} // namespace Tpetra
