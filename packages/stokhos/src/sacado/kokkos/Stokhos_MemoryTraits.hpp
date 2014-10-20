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

#ifndef STOKHOS_MEMORY_TRAITS_HPP
#define STOKHOS_MEMORY_TRAITS_HPP

#include <cstdlib>

#include "Kokkos_Core_fwd.hpp"

// Currently always aligning
#define STOKHOS_ALIGN_MEMORY 1

namespace Stokhos {

//! Traits class encapsulting memory alignment
template <typename MemorySpace>
struct MemoryTraits {

  //! Bytes to which memory allocations are aligned
  static const unsigned Alignment = 8;
  KOKKOS_INLINE_FUNCTION

  //! Allocate aligned memory of given size
  static void* alloc(const size_t size) { return operator new(size); }

  //! Free memory allocated by alloc()
  KOKKOS_INLINE_FUNCTION
  static void free(void *ptr) { operator delete(ptr); }
};

//! Specialization of MemoryTraits for host memory spaces
template <>
struct MemoryTraits< Kokkos::HostSpace > {

//! Bytes to which memory allocations are aligned
#if STOKHOS_ALIGN_MEMORY
#if defined(__MIC__)
  static const unsigned Alignment = 64;
#elif defined(__AVX__)
  static const unsigned Alignment = 32;
#elif defined(__SSE2__)
  static const unsigned Alignment = 16;
#else
  static const unsigned Alignment = 8;
#endif
#else
  static const unsigned Alignment = 8;
#endif

  //! Allocate aligned memory
  /*!
   * Note:  We don't use mm_malloc or posix_memalign, because even though this
   * implementation is host-only, it is potentially callable from host functions
   * marked as device functions (via the KOKKOS_INLINE_FUNCTION maco).
   *
   * Also, we can't call new/delete as we may replace those with a version
   * that calls this.
   */
  KOKKOS_INLINE_FUNCTION
  static void* alloc(const size_t size) {
    void* ptr = 0;
    if (size > 0) {
#if STOKHOS_ALIGN_MEMORY
      const size_t mask = Alignment-1;
      const size_t total_size = size + mask + sizeof(std::ptrdiff_t);
      char *ptr_alloc = reinterpret_cast<char*>(std::malloc(total_size));
      char *ptr_storage = ptr_alloc + sizeof(std::ptrdiff_t);
      char *ptr_body = reinterpret_cast<char*>(
        ( reinterpret_cast<size_t>(ptr_storage) + mask ) & ~mask );
      char *ptr_header = ptr_body - sizeof(std::ptrdiff_t);
      const std::ptrdiff_t offset = ptr_body - ptr_alloc;
      *reinterpret_cast<std::ptrdiff_t*>(ptr_header) = offset;
      ptr = reinterpret_cast<void*>(ptr_body);
#else
      ptr = operator new(size);
#endif
    }
    return ptr;
  }

  //! Free memory allocated by alloc()
  KOKKOS_INLINE_FUNCTION
  static void free(void *ptr) {
    if (ptr != 0) {
#if STOKHOS_ALIGN_MEMORY
      void *ptr_header = reinterpret_cast<char*>(ptr) - sizeof(std::ptrdiff_t);
      const std::ptrdiff_t offset = *reinterpret_cast<std::ptrdiff_t*>(ptr_header);
      void *ptr_alloc = reinterpret_cast<char*>(ptr) - offset;
      std::free(ptr_alloc);
#else
      operator delete(ptr);
#endif
    }
  }
};

//! An aligned STL allocator
template <typename T>
class aligned_allocator {
public:

  typedef T         value_type;
  typedef T*        pointer;
  typedef const T*  const_pointer;
  typedef T&        reference;
  typedef const T&  const_reference;
  typedef size_t    size_type;
  typedef std::ptrdiff_t difference_type;

  typedef Stokhos::MemoryTraits< Kokkos::HostSpace > Traits;

  template <class U> struct rebind { typedef aligned_allocator<U> other; };

  aligned_allocator() {}

  template <class U> aligned_allocator(const aligned_allocator<U>&) {}

  size_type max_size() const {
    return (size_type(~0) - size_type(Traits::Alignment)) / sizeof(T);
  }

  pointer address(reference x) const { return &x; }

  const_pointer address(const_reference x) const { return &x; }

  pointer allocate(size_type n, const void* = 0) {
    size_type count = n * sizeof(T);
    void* ptr = Traits::alloc(count);
    if (ptr == 0) throw std::bad_alloc();
    return reinterpret_cast<pointer>(ptr);
  }

  void deallocate(pointer p, size_type) { Traits::free(p); }

  void construct(pointer p, const_reference val) { new (p) T(val); }

  void destroy(pointer p) { ((T*)p)->~T(); }
};

//! An aligned STL allocator
template <typename T>
class aligned_allocator< const T > {
public:

  typedef T         value_type;
  typedef const T*  pointer;
  typedef const T*  const_pointer;
  typedef const T&  reference;
  typedef const T&  const_reference;
  typedef size_t    size_type;
  typedef std::ptrdiff_t difference_type;

  typedef Stokhos::MemoryTraits< Kokkos::HostSpace > Traits;

  template <class U> struct rebind { typedef aligned_allocator<U> other; };

  aligned_allocator() {}

  template <class U> aligned_allocator(const aligned_allocator<U>&) {}

  size_type max_size() const {
    return (size_type(~0) - size_type(Traits::Alignment)) / sizeof(T);
  }

  const_pointer address(const_reference x) const { return &x; }

  pointer allocate(size_type n, const void* = 0) {
    size_type count = n * sizeof(T);
    void* ptr = Traits::alloc(count);
    if (ptr == 0)  throw std::bad_alloc();
    return reinterpret_cast<pointer>(ptr);
  }

  void deallocate(pointer p, size_type) { Traits::free(p); }

  void construct(pointer p, const_reference val) { new (p) T(val); }

  void destroy(pointer p) { ((T*)p)->~T(); }
};

template <typename T, typename U>
inline bool
operator == (const aligned_allocator<T>&, const aligned_allocator<U>&)
{ return true; }

template <typename T, typename U>
inline bool
operator != (const aligned_allocator<T>&, const aligned_allocator<U>&)
{ return false; }

} // namespace Stokhos

#endif // STOKHOS_MEMORY_TRAITS_HPP
