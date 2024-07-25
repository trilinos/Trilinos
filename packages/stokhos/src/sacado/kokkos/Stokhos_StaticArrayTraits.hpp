// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_STATIC_ARRAY_TRAITS_HPP
#define STOKHOS_STATIC_ARRAY_TRAITS_HPP

#include <cstring>

#include "Sacado_Traits.hpp"

#include "Kokkos_Macros.hpp"

#include "Stokhos_MemoryTraits.hpp"

namespace Stokhos {

  /*!
   * \brief Static array allocation class
   */
  template <typename T, typename device,
            bool isScalar = Sacado::IsScalarType<T>::value>
  struct StaticArrayTraits {};

  /*!
   * \brief Static array allocation class that works for any type
   */
  template <typename T, typename D>
  struct StaticArrayTraits<T, D, false> {

    typedef T value_type;
    typedef D execution_space;

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, volatile T*  dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, T*  dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, volatile T*  dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, T*  dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Zero out array \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void zero(T* dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Zero out array \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void zero(volatile T* dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(T* dest, std::size_t sz, const T& v) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(volatile T* dest, std::size_t sz, const T& v) {
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

  };

  /*!
   * \brief Static array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T, typename D>
  struct StaticArrayTraits<T,D,true> {

    typedef T value_type;
    typedef D execution_space;

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, volatile T* dest, std::size_t sz) {
      // if (sz > 0)
      //   std::memcpy(const_cast<const T*>(dest),const_cast<T*>(src),sz*sizeof(T));
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, T* dest, std::size_t sz) {
      // if (sz > 0)
      //   std::memcpy(dest,const_cast<const T*>(src),sz*sizeof(T));
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, volatile T* dest, std::size_t sz) {
      // if (sz > 0)
      //   std::memcpy(const_cast<T*>(dest),src,sz*sizeof(T));
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, T* dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Zero out array \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void zero(T* dest, std::size_t sz) {
      // if (sz > 0) std::memset(dest,0,sz*sizeof(T));
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Zero out array \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void zero(volatile T* dest, std::size_t sz) {
      // if (sz > 0) std::memset(const_cast<T*>(dest),0,sz*sizeof(T));
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(T* dest, std::size_t sz, T v) {
      //std::memset(dest,v,sz*sizeof(T)); // memset doesn't work if v != 0?
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(volatile T* dest, std::size_t sz, T v) {
      //std::memset(dest,v,sz*sizeof(T)); // memset doesn't work if v != 0?
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

  };

} // namespace Stokhos

#endif // STOKHOS_STATIC_ARRAY_TRAITS_HPP
