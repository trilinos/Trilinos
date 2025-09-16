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
// 

#ifndef STK_IMPL_FIELD_DATA_ALLOCATOR_HPP
#define STK_IMPL_FIELD_DATA_ALLOCATOR_HPP

#include <stk_util/stk_config.h>
#include <stk_util/util/string_utils.hpp>
#include "stk_util/util/ReportHandler.hpp"
#include <stk_util/ngp/NgpSpaces.hpp>
#include "Kokkos_Core.hpp"

#if defined(STK_ASAN_IS_ON) && defined(STK_ASAN_FIELD_ACCESS)
  #include <sanitizer/asan_interface.h>
  constexpr unsigned STK_ASAN_FIELD_PADDING_SIZE = 1024;
#else
  #define ASAN_POISON_MEMORY_REGION(addr, size) ((void)(addr), (void)(size))
  #define ASAN_UNPOISON_MEMORY_REGION(addr, size) ((void)(addr), (void)(size))
  constexpr unsigned STK_ASAN_FIELD_PADDING_SIZE = 0;
#endif

// If a full SIMD access occurs at the end of an irregularly-sized Bucket, we need to
// tolerate accessing uninitialized values as long as they do not walk outside the
// allocation.  Pad the allocation so that it is always a full increment of the maximum
// SIMD size.  This only really has an affect on unusually-small Buckets or Buckets with
// a non-power-of-2 capacity.  This may be smaller than the actual memory alignment size.
//
constexpr unsigned STK_ALIGNMENT_PADDING_SIZE = 64;  // For avx512 SIMD support

//We used to assert that KOKKOS_MEMORY_ALIGNMENT is at least as big as
//STK_ALIGNMENT_PADDING_SIZE.
//That Kokkos macro has been deprecated (as of the 4.7 or so version). Kokkos documentation
//simply indicates that Kokkos Views are 64-byte aligned. There doesn't appear to be a way
//to set Kokkos alignment to a value different than that, nor is there a way to query Kokkos
//View alignment.

namespace stk {

template <typename T>
struct SerialFieldDataAllocator
{
  using DataType = T;
  using HostMemorySpace = stk::ngp::HostMemSpace;
  using DeviceMemorySpace = stk::ngp::MemSpace;
  using HostAllocationType = Kokkos::View<T*, HostMemorySpace>;
  using DeviceAllocationType = Kokkos::View<T*, DeviceMemorySpace>;

  static_assert(sizeof(T) == 1u, "Allocator datatype must be of byte-equivalent length.");

  HostAllocationType host_allocate(size_t size) const {
    const size_t asanSize = size + STK_ASAN_FIELD_PADDING_SIZE;
    auto allocation = HostAllocationType(Kokkos::view_alloc("HostFieldData", Kokkos::WithoutInitializing), asanSize);
    ASAN_POISON_MEMORY_REGION(allocation.data(), asanSize);
    return allocation;
  }

  DeviceAllocationType device_allocate(size_t size) const {
    auto allocation = DeviceAllocationType(Kokkos::view_alloc("DeviceFieldData", Kokkos::WithoutInitializing), size);
    return allocation;
  }

  template <typename PtrType>
  PtrType* get_host_pointer_for_device(PtrType* hostPointer) const {
    return hostPointer;
  }
};

#if defined(KOKKOS_ENABLE_CUDA)
template <typename T>
struct CudaFieldDataAllocator
{
  using DataType = T;
  using HostMemorySpace = stk::ngp::HostPinnedSpace;
  using DeviceMemorySpace = stk::ngp::MemSpace;
  using HostAllocationType = Kokkos::View<T*, HostMemorySpace>;
  using DeviceAllocationType = Kokkos::View<T*, DeviceMemorySpace>;

  static_assert(sizeof(T) == 1u, "Allocator datatype must be of byte-equivalent length.");

  HostAllocationType host_allocate(size_t size) const {
    return HostAllocationType(Kokkos::view_alloc("HostFieldData", Kokkos::WithoutInitializing), size);
  }

  DeviceAllocationType device_allocate(size_t size) const {
    return DeviceAllocationType(Kokkos::view_alloc("DeviceFieldData", Kokkos::WithoutInitializing), size);
  }

  template <typename PtrType>
  PtrType* get_host_pointer_for_device(PtrType* hostPointer) const {
    PtrType* translatedHostPointer = hostPointer;
    if (hostPointer != nullptr) {
      cudaError_t status = cudaHostGetDevicePointer((void**)&translatedHostPointer, (void*)hostPointer, 0);
      STK_ThrowRequireMsg(status == cudaSuccess, "Something went wrong during cudaHostGetDevicePointer(): " +
                          std::string(cudaGetErrorString(status)));
    }
    return translatedHostPointer;
  }
};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <typename T>
class HipFieldDataAllocator
{
public:
  using DataType = T;
#if defined(STK_UNIFIED_MEMORY)
  using HostMemorySpace = std::conditional_t<stk::ngp::DeviceAccessibleFromHost,
                                             stk::ngp::MemSpace,          // Unified memory, so allocate on device (for performance)
                                             stk::ngp::HostPinnedSpace>;  // Non-unified memory, so allocate on host
#else
  using HostMemorySpace = stk::ngp::HostPinnedSpace;
#endif
  using DeviceMemorySpace = stk::ngp::MemSpace;
  using HostAllocationType = Kokkos::View<T*, HostMemorySpace>;
  using DeviceAllocationType = Kokkos::View<T*, DeviceMemorySpace>;

  static_assert(sizeof(T) == 1u, "Allocator datatype must be of byte-equivalent length.");

  HostAllocationType host_allocate(size_t size) const {
    return HostAllocationType(Kokkos::view_alloc("HostFieldData", Kokkos::WithoutInitializing), size);
  }

  DeviceAllocationType device_allocate(size_t size) const {
    return DeviceAllocationType(Kokkos::view_alloc("DeviceFieldData", Kokkos::WithoutInitializing), size);
  }

  template <typename PtrType>
  PtrType* get_host_pointer_for_device(PtrType* hostPointer) const {
    PtrType* translatedHostPointer = hostPointer;
    if (hostPointer != nullptr) {
#if defined(STK_UNIFIED_MEMORY)
      if constexpr (not stk::ngp::DeviceAccessibleFromHost) {
        hipError_t status = hipHostGetDevicePointer((void**)&translatedHostPointer, (void*)hostPointer, 0);
        STK_ThrowRequireMsg(status == hipSuccess, "Something went wrong during hipHostGetDevicePointer(): " +
                            std::string(hipGetErrorString(status)));
      }
#else
      hipError_t status = hipHostGetDevicePointer((void**)&translatedHostPointer, (void*)hostPointer, 0);
      STK_ThrowRequireMsg(status == hipSuccess, "Something went wrong during hipHostGetDevicePointer(): " +
                          std::string(hipGetErrorString(status)));
#endif
    }
    return translatedHostPointer;
  }
};
#endif


#if defined(KOKKOS_ENABLE_CUDA)

template <typename DataType>
using FieldDataAllocator = CudaFieldDataAllocator<DataType>;

#elif defined(KOKKOS_ENABLE_HIP)

template <typename DataType>
using FieldDataAllocator = HipFieldDataAllocator<DataType>;

#else

template <typename DataType>
using FieldDataAllocator = SerialFieldDataAllocator<DataType>;

#endif

}  // namespace stk

#endif
