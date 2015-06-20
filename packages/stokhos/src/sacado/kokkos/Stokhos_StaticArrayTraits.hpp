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

#ifndef STOKHOS_STATIC_ARRAY_TRAITS_HPP
#define STOKHOS_STATIC_ARRAY_TRAITS_HPP

#include <cstring>

#include "Sacado_Traits.hpp"

#include "Kokkos_Macros.hpp"

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
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, T* dest, std::size_t sz) {
      // if (sz > 0)
      //   std::memcpy(dest,const_cast<const T*>(src),sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, volatile T* dest, std::size_t sz) {
      // if (sz > 0)
      //   std::memcpy(const_cast<T*>(dest),src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, T* dest, std::size_t sz) {
      if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
    }

    //! Zero out array \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void zero(T* dest, std::size_t sz) {
      if (sz > 0) std::memset(dest,0,sz*sizeof(T));
    }

    //! Zero out array \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void zero(volatile T* dest, std::size_t sz) {
      // if (sz > 0) std::memset(const_cast<T*>(dest),0,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(T* dest, std::size_t sz, T v) {
      //std::memset(dest,v,sz*sizeof(T)); // memset doesn't work if v != 0?
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(volatile T* dest, std::size_t sz, T v) {
      //std::memset(dest,v,sz*sizeof(T)); // memset doesn't work if v != 0?
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

  };

} // namespace Stokhos

#endif // STOKHOS_STATIC_ARRAY_TRAITS_HPP
