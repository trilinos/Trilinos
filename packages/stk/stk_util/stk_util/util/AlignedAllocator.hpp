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

#ifndef STK_ALIGNED_ALLOCATOR_HPP
#define STK_ALIGNED_ALLOCATOR_HPP

#include "stk_util/util/AdjustForAlignment.hpp"

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include <cstdlib>
#include <stdlib.h>
#include <memory>

#include <Kokkos_Macros.hpp>

namespace non_std {

template <typename T, size_t Alignment>
class BaseAlignedAllocator
{
public:
  using ValueType = T;
  using Pointer = T*;
  using ConstPointer = const T*;
  using SizeType = std::size_t;

  template <typename U, size_t A>
  struct rebind
  {
    using OtherType = BaseAlignedAllocator<U,A>;
  };

  BaseAlignedAllocator() {}

  BaseAlignedAllocator(const BaseAlignedAllocator&) {}

  template <typename U, size_t A>
  BaseAlignedAllocator (const BaseAlignedAllocator <U,A>&) {}

  ~BaseAlignedAllocator() {}
};

template <class T, size_t Alignment>
struct HostAlignedAllocator
  : public std::allocator<T>, public BaseAlignedAllocator<T, Alignment>
{
  using Pointer = typename BaseAlignedAllocator<T, Alignment>::Pointer;
  using ConstPointer = typename BaseAlignedAllocator<T, Alignment>::ConstPointer;
  using SizeType = typename BaseAlignedAllocator<T, Alignment>::SizeType;

  template <class U>
  struct rebind { typedef HostAlignedAllocator<U, Alignment> other; };

  HostAlignedAllocator() noexcept
  {
    static_assert((Alignment & (Alignment - 1UL)) == 0, "Alignment must be power of 2");
  }

  HostAlignedAllocator(const HostAlignedAllocator& other) noexcept
      : std::allocator<T>(other), BaseAlignedAllocator<T, Alignment>(other)
  {
  }

  template <class U>
  HostAlignedAllocator(const HostAlignedAllocator<U,Alignment>&) noexcept { }

  ~HostAlignedAllocator() noexcept { }

  inline Pointer allocate(SizeType n) {
    return allocate(n, ConstPointer(0));
  }
    
  inline Pointer allocate(SizeType n, ConstPointer) {
    size_t requested = n * sizeof(T);
    void* p = std::aligned_alloc(Alignment, stk::adjust_up_to_alignment_boundary(requested, Alignment));
    if (!p) {
      p = nullptr;
      throw std::bad_alloc();
    }
    return static_cast<Pointer>(p);
  }

  inline void deallocate(Pointer p, SizeType) {
    free(p);
  }
};

#ifdef KOKKOS_ENABLE_CUDA
namespace impl {
  void* CUDAPinnedAndMappedAlignedAllocate(size_t size);
  void CUDAPinnedAndMappedAlignedDeallocate(void* p);
}

template <typename T, size_t Alignment>
class CUDAPinnedAndMappedAlignedAllocator
  : public std::allocator<T>, public BaseAlignedAllocator<T, Alignment>
{
public:
  using ValueType = typename BaseAlignedAllocator<T, Alignment>::ValueType;
  using Pointer = typename BaseAlignedAllocator<T, Alignment>::Pointer;
  using SizeType = typename BaseAlignedAllocator<T, Alignment>::SizeType;

  template <class U>
  struct rebind { typedef CUDAPinnedAndMappedAlignedAllocator<U, Alignment> other; };

  CUDAPinnedAndMappedAlignedAllocator() {}

  CUDAPinnedAndMappedAlignedAllocator(const CUDAPinnedAndMappedAlignedAllocator&) {}

  template <typename U>
  CUDAPinnedAndMappedAlignedAllocator (const CUDAPinnedAndMappedAlignedAllocator<U, Alignment>&) {}

  ~CUDAPinnedAndMappedAlignedAllocator() {}

  inline Pointer allocate(SizeType num, const void* = 0)
  {
    size_t size = num * sizeof(ValueType);

    return reinterpret_cast<T*>(impl::CUDAPinnedAndMappedAlignedAllocate(size));
  }

  inline void deallocate(Pointer p, SizeType)
  {
    impl::CUDAPinnedAndMappedAlignedDeallocate(p);
  }
};
#endif

#ifdef KOKKOS_ENABLE_HIP
namespace impl {
  void* HIPPinnedAndMappedAlignedAllocate(size_t size);
  void HIPPinnedAndMappedAlignedDeallocate(void* p);
}

template <typename T, size_t Alignment>
class HIPPinnedAndMappedAlignedAllocator
  : public std::allocator<T>, public BaseAlignedAllocator<T, Alignment>
{
public:
  using ValueType = typename BaseAlignedAllocator<T, Alignment>::ValueType;
  using Pointer = typename BaseAlignedAllocator<T, Alignment>::Pointer;
  using SizeType = typename BaseAlignedAllocator<T, Alignment>::SizeType;

  template <class U>
  struct rebind { typedef HIPPinnedAndMappedAlignedAllocator<U, Alignment> other; };

  HIPPinnedAndMappedAlignedAllocator() {}

  HIPPinnedAndMappedAlignedAllocator(const HIPPinnedAndMappedAlignedAllocator&) {}

  template <typename U>
  HIPPinnedAndMappedAlignedAllocator (const HIPPinnedAndMappedAlignedAllocator<U, Alignment>&) {}

  ~HIPPinnedAndMappedAlignedAllocator() {}

  inline Pointer allocate(SizeType num, const void* = 0)
  {
    size_t size = num * sizeof(ValueType);

    return reinterpret_cast<T*>(impl::HIPPinnedAndMappedAlignedAllocate(size));
  }

  inline void deallocate(Pointer p, SizeType)
  {
    impl::HIPPinnedAndMappedAlignedDeallocate(p);
  }
};
#endif

template <class T1, size_t A1, class T2, size_t A2> inline
bool operator == (const BaseAlignedAllocator<T1,A1> &, const BaseAlignedAllocator<T2,A2> &) {
  return true;
}

template <class T1, size_t A1, class T2, size_t A2> inline
bool operator != (const BaseAlignedAllocator<T1,A1> &, const BaseAlignedAllocator<T2,A2> &) {
  return false;
}

#ifdef KOKKOS_ENABLE_CUDA
template <typename ValueType, size_t Alignment>
using AlignedAllocator = CUDAPinnedAndMappedAlignedAllocator<ValueType, Alignment>;
#elif defined(KOKKOS_ENABLE_HIP)
template <typename ValueType, size_t Alignment>
using AlignedAllocator = HIPPinnedAndMappedAlignedAllocator<ValueType, Alignment>;
#else
template <typename ValueType, size_t Alignment>
using AlignedAllocator = HostAlignedAllocator<ValueType, Alignment>;
#endif

} // namespace non_std

#endif
