// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MEMORY_TRAITS_HPP
#define STOKHOS_MEMORY_TRAITS_HPP

#include <cstdlib>

#include "Kokkos_Core_fwd.hpp"

// Currently always aligning
#define STOKHOS_ALIGN_MEMORY 1

// Uncomment this if you know all accesses will be aligned.  Tthis is true
// if all Stokhos variables are coming from Kokkos allocations, new/delete
// or stack variables.  However it may not be true for C allocations (e.g.,
// through MPI).
// #define STOKHOS_ASSUME_ALIGNED

// ivdep is necessary to get the intel compiler to vectorize through
// expresion template assignent operators
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#define STOKHOS_HAVE_PRAGMA_IVDEP
#endif

// unrolling appears to slow everything down
#if 0 && ( defined(__INTEL_COMPILER) || defined(__CUDA_ARCH__) )
#define STOKHOS_HAVE_PRAGMA_UNROLL
#endif

// assume all memory accesses are aligned appropriately for aligned
// vector loads
#if defined(STOKHOS_ALIGN_MEMORY) && defined(STOKHOS_ASSUME_ALIGNED) && defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#define STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#endif

namespace Stokhos {

//! Traits class encapsulting memory alignment
template <typename MemorySpace>
struct MemoryTraits {

  //! Bytes to which memory allocations are aligned
  static const unsigned Alignment = 8;

  //! Allocate aligned memory of given size
  KOKKOS_INLINE_FUNCTION
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
