// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_DYNAMICARRAYTRAITS_HPP
#define SACADO_DYNAMICARRAYTRAITS_HPP

#include <new>
#include <cstring>

#include "Sacado_Traits.hpp"
#if defined(HAVE_SACADO_KOKKOSCORE)
#include "Kokkos_Core.hpp"
#endif

namespace Sacado {

  /*!
   * \brief Dynamic array allocation class that works for any type
   */
  template <typename T, bool isScalar = IsScalarType<T>::value>
  struct ds_array {

    KOKKOS_INLINE_FUNCTION
    static T* my_alloc(const int sz) {
#if defined(__CUDACC__) && defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
      T* m;
      cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
      return m;
    }

    KOKKOS_INLINE_FUNCTION
    static void my_free(T* m, int sz) {
      if (sz > 0) {
#if defined(__CUDACC__) && defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaFree(m);
#else
        operator delete((void*) m);
#endif
      }
    }

    //! Get memory for new array of length \c sz
    KOKKOS_INLINE_FUNCTION
    static T* get(int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        T* p = m;
        for (int i=0; i<sz; ++i)
          new (p++) T();
        return m;
      }
      return NULL;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    KOKKOS_INLINE_FUNCTION
    static T* get_and_fill(int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        T* p = m;
        for (int i=0; i<sz; ++i)
          new (p++) T(0.0);
        return m;
      }
      return NULL;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    KOKKOS_INLINE_FUNCTION
    static T* get_and_fill(const T* src, int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        T* p = m;
        for (int i=0; i<sz; ++i)
          new (p++) T(*(src++));
        return m;
      }
      return NULL;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    KOKKOS_INLINE_FUNCTION
    static T* strided_get_and_fill(const T* src, int stride, int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        T* p = m;
        for (int i=0; i<sz; ++i) {
          new (p++) T(*(src));
          src += stride;
        }
        return m;
      }
      return NULL;
    }

    //! Copy array from \c src to \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void copy(const T* src, T*  dest, int sz) {
      for (int i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void strided_copy(const T* src, int src_stride,
                                    T* dest, int dest_stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = *(src);
        dest += dest_stride;
        src += src_stride;
      }
    }

    //! Zero out array \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
      for (int i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Zero out array \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void strided_zero(T* dest, int stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = T(0.);
        dest += stride;
      }
    }

    //! Destroy array elements and release memory
    KOKKOS_INLINE_FUNCTION
    static void destroy_and_release(T* m, int sz) {
      T* e = m+sz;
      for (T* b = m; b!=e; b++)
        b->~T();
      my_free(m, sz);
    }
  };

  /*!
   * \brief Dynamic array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T>
  struct ds_array<T,true> {

    KOKKOS_INLINE_FUNCTION
    static T* my_alloc(const int sz) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
      T* m;
      cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
      return m;
    }

    KOKKOS_INLINE_FUNCTION
    static void my_free(T* m, int sz) {
      if (sz > 0) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaFree(m);
#else
        operator delete((void*) m);
#endif
      }
    }

    //! Get memory for new array of length \c sz
    KOKKOS_INLINE_FUNCTION
    static T* get(int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        return m;
      }
      return NULL;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    KOKKOS_INLINE_FUNCTION
    static T* get_and_fill(int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        std::memset(m,0,sz*sizeof(T));
        return m;
      }
      return NULL;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    KOKKOS_INLINE_FUNCTION
    static T* get_and_fill(const T* src, int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        for (int i=0; i<sz; ++i)
          m[i] = src[i];
        return m;
      }
      return NULL;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    KOKKOS_INLINE_FUNCTION
    static T* strided_get_and_fill(const T* src, int stride, int sz) {
      if (sz > 0) {
        T* m = my_alloc(sz);
        for (int i=0; i<sz; ++i)
          m[i] = src[i*stride];
        return m;
      }
      return NULL;
    }

    //! Copy array from \c src to \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void copy(const T* src, T* dest, int sz) {
      std::memcpy(dest,src,sz*sizeof(T));
    }

    //! Copy array from \c src to \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void strided_copy(const T* src, int src_stride,
                                    T* dest, int dest_stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = *(src);
        dest += dest_stride;
        src += src_stride;
      }
    }

    //! Zero out array \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
      if (sz > 0)
        std::memset(dest,0,sz*sizeof(T));
    }

    //! Zero out array \c dest of length \c sz
    KOKKOS_INLINE_FUNCTION
    static void strided_zero(T* dest, int stride, int sz) {
      for (int i=0; i<sz; ++i) {
        *(dest) = T(0.);
        dest += stride;
      }
    }

    //! Destroy array elements and release memory
    KOKKOS_INLINE_FUNCTION
    static void destroy_and_release(T* m, int sz) {
      my_free(m, sz);
    }
  };

} // namespace Sacado

#endif // SACADO_DYNAMICARRAY_HPP
