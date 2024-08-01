// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#ifdef KOKKOS_ENABLE_HIP

void* hip_memory_ = nullptr;
size_t hip_memory_size_ = 0;

void finalize_hip_memory ()
{
  if (hip_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::HIPSpace> (hip_memory_);
    hip_memory_ = nullptr;
    hip_memory_size_ = 0;
  }
}

void* hip_host_pinned_memory_ = nullptr;
size_t hip_host_pinned_memory_size_ = 0;

void finalize_hip_host_pinned_memory ()
{
  if (hip_host_pinned_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::HIPHostPinnedSpace> (hip_host_pinned_memory_);
    hip_host_pinned_memory_ = nullptr;
    hip_host_pinned_memory_size_ = 0;
  }
}
#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL

void* sycl_memory_ = nullptr;
size_t sycl_memory_size_ = 0;

void finalize_sycl_memory ()
{
  if (sycl_memory_ != nullptr) {
    Kokkos::kokkos_free<Kokkos::Experimental::SYCLDeviceUSMSpace> (sycl_memory_);
    sycl_memory_ = nullptr;
    sycl_memory_size_ = 0;
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

#ifdef KOKKOS_ENABLE_HIP

void*
StaticKokkosAllocation<Kokkos::HIPSpace>::
resize (Kokkos::HIPSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::HIPSpace;
  static bool created_finalize_hook = false;

  if (size > hip_memory_size_) {
    if (hip_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (hip_memory_);
    }
    const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
    hip_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    hip_memory_size_ = size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_hip_memory);
    created_finalize_hook = true;
  }

  return hip_memory_;
}

void*
StaticKokkosAllocation<Kokkos::HIPHostPinnedSpace>::
resize (Kokkos::HIPHostPinnedSpace /* space */,
        const size_t size)
{
  using memory_space = Kokkos::HIPHostPinnedSpace;
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
  static bool created_finalize_hook = false;

  if (size > sycl_memory_size_) {
    if (sycl_memory_ != nullptr) {
      Kokkos::kokkos_free<memory_space> (sycl_memory_);
    }
    const size_t req_size = size > minimum_initial_size ? size : minimum_initial_size;
    sycl_memory_ = Kokkos::kokkos_malloc<memory_space> (req_size);
    sycl_memory_size_ = size;
  }
  if (! created_finalize_hook) {
    Kokkos::push_finalize_hook (finalize_sycl_memory);
    created_finalize_hook = true;
  }

  return sycl_memory_;
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
