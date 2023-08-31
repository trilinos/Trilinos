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
bool created_cuda_finalize_hook = false;
Kokkos::Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space> * cuda_pool_=nullptr;
void* cuda_memory_ = nullptr;
size_t cuda_memory_size_ = 0;

void finalize_cuda_memory ()
{
  if (cuda_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::CudaSpace> (cuda_memory_);
    cuda_memory_ = nullptr;
    cuda_memory_size_ = 0;
  }

  if(cuda_pool_ != nullptr) {
    delete cuda_pool_;
    cuda_pool_ = nullptr;
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

#ifdef KOKKOS_ENABLE_HIP
bool created_hip_finalize_hook = false;
Kokkos::Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space> * hip_pool_=nullptr;
void* hip_memory_ = nullptr;
size_t hip_memory_size_ = 0;

void finalize_hip_memory ()
{
  if (hip_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::Experimental::HIPSpace> (hip_memory_);
    hip_memory_ = nullptr;
    hip_memory_size_ = 0;
  }

  if(hip_pool_ != nullptr) {
    delete hip_pool_;
    hip_pool_ = nullptr;
  }    
}

void* hip_host_pinned_memory_ = nullptr;
size_t hip_host_pinned_memory_size_ = 0;

void finalize_hip_host_pinned_memory ()
{
  if (hip_host_pinned_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::Experimental::HIPHostPinnedSpace> (hip_host_pinned_memory_);
    hip_host_pinned_memory_ = nullptr;
    hip_host_pinned_memory_size_ = 0;
  }
}
#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL
bool created_sycl_finalize_hook = false;
Kokkos::Random_XorShift64_Pool<typename Kokkos::SYCLDeviceUSMSpace::execution_space> * sycl_pool_=nullptr;
void* sycl_memory_ = nullptr;
size_t sycl_memory_size_ = 0;

void finalize_sycl_memory ()
{
  if (sycl_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::Experimental::SYCLDeviceUSMSpace> (sycl_memory_);
    sycl_memory_ = nullptr;
    sycl_memory_size_ = 0;
  }

  if(sycl_pool_ != nullptr) {
    delete sycl_pool_;
    sycl_pool_ = nullptr;
  }    
}

void* sycl_shared_memory_ = nullptr;
size_t sycl_shared_memory_size_ = 0;

void finalize_sycl_shared_memory ()
{
  if (sycl_shared_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::Experimental::SYCLSharedUSMSpace> (sycl_shared_memory_);
    sycl_shared_memory_ = nullptr;
    sycl_shared_memory_size_ = 0;
  }
}
#endif // KOKKOS_ENABLE_SYCL

bool created_host_finalize_hook = false;
Kokkos::Random_XorShift64_Pool<typename Kokkos::HostSpace::execution_space> * host_pool_=nullptr;

void* host_memory_ = nullptr;
size_t host_memory_size_ = 0;

void finalize_host_memory ()
{
  if (host_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::HostSpace> (host_memory_);
    host_memory_ = nullptr;
    host_memory_size_ = 0;
  }

  if(host_pool_ != nullptr) {
    delete host_pool_;
    host_pool_ = nullptr;
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
  if (size > cuda_memory_size_) {
    if (cuda_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (cuda_memory_);
    }
    const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
    cuda_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    cuda_memory_size_ = size;
  }
  if (! created_cuda_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_cuda_memory);
    created_cuda_finalize_hook = true;
  }

  return cuda_memory_;
}

Kokkos::Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space> &
Static_Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space>::
getPool(unsigned int seed) {
  using pool_type = Kokkos::Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space>;

  if(cuda_pool_ == nullptr) {
    cuda_pool_ = new pool_type(seed);
    if (! created_cuda_finalize_hook) {
      Kokkos::push_finalize_hook (finalize_cuda_memory);
      created_cuda_finalize_hook = true;
    }
  }
  return *cuda_pool_;
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

#ifdef KOKKOS_ENABLE_HIP

void*
StaticKokkosAllocation<Kokkos::Experimental::HIPSpace>::
resize (Kokkos::Experimental::HIPSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::Experimental::HIPSpace;

  if (size > hip_memory_size_) {
    if (hip_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (hip_memory_);
    }
    const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
    hip_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    hip_memory_size_ = size;
  }
  if (! created_hip_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_hip_memory);
    created_hip_finalize_hook = true;
  }

  return hip_memory_;
}

Kokkos::Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space> &
Static_Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space>::
getPool(unsigned int seed) {
  using pool_type = Kokkos::Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space>;

  if(hip_pool_ == nullptr) {
    hip_pool_ = new pool_type(seed);
    if (! created_hip_finalize_hook) {
      Kokkos::push_finalize_hook (finalize_cuda_memory);
      created_hip_finalize_hook = true;
    }
  }
  return *hip_pool_;
}

void*
StaticKokkosAllocation<Kokkos::Experimental::HIPHostPinnedSpace>::
resize (Kokkos::Experimental::HIPHostPinnedSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::Experimental::HIPHostPinnedSpace;
  static bool created_finalize_hook = false;

  const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
  if (req_size > hip_host_pinned_memory_size_) {
    if (hip_host_pinned_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (hip_host_pinned_memory_);
    }
    hip_host_pinned_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    hip_host_pinned_memory_size_ = req_size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_hip_host_pinned_memory);
    created_finalize_hook = true;
  }

  return hip_host_pinned_memory_;
}

#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL

template <>
void*
StaticKokkosAllocation<Kokkos::Experimental::SYCLDeviceUSMSpace>::
resize (Kokkos::Experimental::SYCLDeviceUSMSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::Experimental::SYCLDeviceUSMSpace;

  if (size > sycl_memory_size_) {
    if (sycl_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (sycl_memory_);
    }
    const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
    sycl_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    sycl_memory_size_ = size;
  }
  if (! created_syck_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_sycl_memory);
    created_sycl)finalize_hook = true;
  }

  return sycl_memory_;
}


Kokkos::Random_XorShift64_Pool<typename Kokkos::SYCLDeviceUSMSpace::execution_space> &
Static_Random_XorShift64_Pool<typename Kokkos::SYCLDeviceUSMSpace::execution_space>::
getPool(unsigned int seed) {
  using pool_type = Kokkos::Random_XorShift64_Pool<typename Kokkos::SYCLDeviceUSMSpace::execution_space>;

  if(sycl_pool_ == nullptr) {
    sycl_pool_ = new pool_type(seed);
    if (! created_sycl_finalize_hook) {
      Kokkos::push_finalize_hook (finalize_cuda_memory);
      created_sycl_finalize_hook = true;
    }
  }
  return *sycl_pool_;
}


template <>
void*
StaticKokkosAllocation<Kokkos::Experimental::SYCLSharedUSMSpace>::
resize (Kokkos::Experimental::SYCLSharedUSMSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::Experimental::SYCLSharedUSMSpace;
  static bool created_finalize_hook = false;

  const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
  if (req_size > sycl_shared_memory_size_) {
    if (sycl_shared_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (sycl_shared_memory_);
    }
    sycl_shared_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    sycl_shared_memory_size_ = req_size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_sycl_shared_memory);
    created_finalize_hook = true;
  }

  return sycl_shared_memory_;
}

#endif // KOKKOS_ENABLE_SYCL

void*
StaticKokkosAllocation<Kokkos::HostSpace>::
resize (Kokkos::HostSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::HostSpace;

  const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
  if (req_size > host_memory_size_) {
    if (host_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (host_memory_);
    }
    host_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    host_memory_size_ = req_size;
  }
  if (! created_host_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_host_memory);
    created_host_finalize_hook = true;
  }

  return host_memory_;
}

Kokkos::Random_XorShift64_Pool<typename Kokkos::HostSpace::execution_space> &
Static_Random_XorShift64_Pool<typename Kokkos::HostSpace::execution_space>::
getPool(unsigned int seed) {
  using pool_type = Kokkos::Random_XorShift64_Pool<typename Kokkos::HostSpace::execution_space>;

  if(host_pool_ == nullptr) {
    host_pool_ = new pool_type(seed);
    if (! created_host_finalize_hook) {
      Kokkos::push_finalize_hook (finalize_cuda_memory);
      created_host_finalize_hook = true;
    }
  }
  return *host_pool_;
}



} // namespace Impl
} // namespace Details
} // namespace Tpetra
