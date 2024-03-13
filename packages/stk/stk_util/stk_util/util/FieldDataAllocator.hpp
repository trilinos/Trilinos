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

#include <cstddef>      // for size_t, ptrdiff_t
#include <cstdlib>      // for malloc, free
#include <limits>       // for numeric_limits
#include <type_traits>  // for is_same
#include <unistd.h>    // for sysconf, _SC_PAGE_SIZE
#include <new>         // for bad_alloc
#include <iostream>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/stk_config.h>
#include "Kokkos_Core.hpp"

#if defined(STK_ASAN_IS_ON) && defined(STK_ASAN_FIELD_ACCESS)
  #include <sanitizer/asan_interface.h>
  constexpr unsigned ASAN_FIELD_PADDING_SIZE = 1024;
#else
  #define ASAN_POISON_MEMORY_REGION(addr, size) ((void)(addr), (void)(size))
  #define ASAN_UNPOISON_MEMORY_REGION(addr, size) ((void)(addr), (void)(size))
  constexpr unsigned ASAN_FIELD_PADDING_SIZE = 0;
#endif

namespace stk {

namespace impl {

constexpr unsigned DEFAULT_FIELD_ALIGNMENT_BYTES = 64;  // For avx512 SIMD support

template <typename T>
class BaseFieldDataAllocator
{
public:
  using ValueType = T;
  using Pointer = T*;
  using ConstPointer = const T*;
  using Reference = T&;
  using ConstReference = const T&;
  using SizeType = std::size_t;
  using DifferenceType = std::ptrdiff_t;

  template <typename U>
  struct rebind
  {
    using OtherType = BaseFieldDataAllocator<U>;
  };

  BaseFieldDataAllocator() {}

  BaseFieldDataAllocator(const BaseFieldDataAllocator&) {}

  template <typename U>
  BaseFieldDataAllocator (const BaseFieldDataAllocator <U>&) {}

  ~BaseFieldDataAllocator() {}

  static Pointer address(Reference value) { return &value; }
  static ConstPointer address(ConstReference value) { return &value; }

  static SizeType max_size()
  {
    return std::numeric_limits<std::size_t>::max() / sizeof(ValueType);
  }

  static void construct(Pointer p, const ValueType& value)
  {
    new(p)ValueType(value);
  }

  static void destroy(Pointer p)
  {
    p->~ValueType();
  }
};

template <typename T, unsigned ALIGNMENT = DEFAULT_FIELD_ALIGNMENT_BYTES>
class AlignedFieldDataAllocator : public BaseFieldDataAllocator<T>
{
public:
  using ValueType = typename BaseFieldDataAllocator<T>::ValueType;
  using Pointer = typename BaseFieldDataAllocator<T>::Pointer;
  using SizeType = typename BaseFieldDataAllocator<T>::SizeType;

  AlignedFieldDataAllocator() = default;

  AlignedFieldDataAllocator(const AlignedFieldDataAllocator&) {}

  template <typename U>
  AlignedFieldDataAllocator (const AlignedFieldDataAllocator<U>&) {}

  ~AlignedFieldDataAllocator() = default;

  static Pointer allocate(SizeType num, const void* = 0)
  {
    size_t alignedSize = ((num * sizeof(ValueType) + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT + ASAN_FIELD_PADDING_SIZE;
    void * allocatedMemory = aligned_alloc(ALIGNMENT, alignedSize);
    ASAN_POISON_MEMORY_REGION(allocatedMemory, alignedSize);
    return static_cast<Pointer>(allocatedMemory);
  }

  static void deallocate(Pointer p, SizeType)
  {
    free(p);
  }
};

#ifdef KOKKOS_ENABLE_CUDA
template <typename T>
class CUDAPinnedAndMappedAllocator : public BaseFieldDataAllocator<T>
{
public:
  using ValueType = typename BaseFieldDataAllocator<T>::ValueType;
  using Pointer = typename BaseFieldDataAllocator<T>::Pointer;
  using SizeType = typename BaseFieldDataAllocator<T>::SizeType;

  CUDAPinnedAndMappedAllocator() {}

  CUDAPinnedAndMappedAllocator(const CUDAPinnedAndMappedAllocator&) {}

  template <typename U>
  CUDAPinnedAndMappedAllocator (const CUDAPinnedAndMappedAllocator<U>&) {}

  ~CUDAPinnedAndMappedAllocator() {}

  static Pointer allocate(SizeType num, const void* = 0)
  {
    size_t size = num * sizeof(ValueType);

    void* ret;
    cudaError_t status = cudaHostAlloc(&ret, size, cudaHostAllocMapped);

    STK_ThrowRequireMsg(status == cudaSuccess, "Error during CUDAPinnedAndMappedAllocator::allocate: " + std::string(cudaGetErrorString(status)));

    return reinterpret_cast<T*>(ret);
  }

  static void deallocate(Pointer p, SizeType)
  {
    cudaFreeHost(p);
  }
};
#endif

#ifdef KOKKOS_ENABLE_HIP
template <typename T>
class HIPPinnedAndMappedAllocator : public BaseFieldDataAllocator<T>
{
public:
  using ValueType = typename BaseFieldDataAllocator<T>::ValueType;
  using Pointer = typename BaseFieldDataAllocator<T>::Pointer;
  using SizeType = typename BaseFieldDataAllocator<T>::SizeType;

  HIPPinnedAndMappedAllocator() {}

  HIPPinnedAndMappedAllocator(const HIPPinnedAndMappedAllocator&) {}

  template <typename U>
  HIPPinnedAndMappedAllocator (const HIPPinnedAndMappedAllocator<U>&) {}

  ~HIPPinnedAndMappedAllocator() {}

  static Pointer allocate(SizeType num, const void* = 0)
  {
    size_t size = num * sizeof(ValueType);

    void* ret;
    hipError_t status = hipHostMalloc(&ret, size, hipHostMallocMapped);

    STK_ThrowRequireMsg(status == hipSuccess, "Error during HIPPinnedAndMappedAllocator::allocate: " + std::string(hipGetErrorString(status)));

    return reinterpret_cast<T*>(ret);
  }

  static void deallocate(Pointer p, SizeType)
  {
    hipError_t status = hipHostFree(p);
    STK_ThrowRequireMsg(status == hipSuccess, "Error during HIPPinnedAndMappedAllocator::deallocate: hipHostFree returned hipError_t=="<<status);
  }

};
#endif

template <typename T1, typename T2>
inline bool operator==(const BaseFieldDataAllocator<T1>&, const BaseFieldDataAllocator<T2>&)
{ return std::is_same<T1,T2>::value; }

template <typename T1, typename T2>
inline bool operator!=(const BaseFieldDataAllocator<T1>&, const BaseFieldDataAllocator<T2>&)
{ return !std::is_same<T1,T2>::value; }

#ifdef KOKKOS_ENABLE_CUDA
template <typename ValueType>
using FieldDataAllocator = CUDAPinnedAndMappedAllocator<ValueType>;
#elif defined(KOKKOS_ENABLE_HIP)
template <typename ValueType>
using FieldDataAllocator = HIPPinnedAndMappedAllocator<ValueType>;
#else
template <typename ValueType>
using FieldDataAllocator = AlignedFieldDataAllocator<ValueType>;
#endif

} // namespace stk
}

#endif
