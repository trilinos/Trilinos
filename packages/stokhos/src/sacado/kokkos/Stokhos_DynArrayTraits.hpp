// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_DYN_ARRAY_TRAITS_HPP
#define STOKHOS_DYN_ARRAY_TRAITS_HPP

#include <new>
#include <cstring>

#include "Kokkos_Core_fwd.hpp"

namespace Stokhos {

  //! Base template specification for %IsScalarType
  /*!
   * The %IsScalarType classes provide a mechanism for computing the
   * determining whether a type is a scalar type (float, double, etc...)
   */
  template <typename T> struct IsScalarType2 {
    static const bool value = false;
  };

  //! Specialization of above classes to built-in types
#define STOKHOS_BUILTIN_SPECIALIZATION(t)                  \
  template <> struct IsScalarType2< t > {                 \
    static const bool value = true;                       \
  };

  STOKHOS_BUILTIN_SPECIALIZATION(float)
  STOKHOS_BUILTIN_SPECIALIZATION(double)
  STOKHOS_BUILTIN_SPECIALIZATION(int)
  STOKHOS_BUILTIN_SPECIALIZATION(long)

#undef STOKHOS_BUILTIN_SPECIALIZATION

  /*!
   * \brief Dynamic array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T, typename device_t,
            bool isScalar = IsScalarType2<T>::value>
  struct DynArrayTraits {

    typedef T value_type;
    typedef device_t execution_space;

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, volatile T*  dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, T*  dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, volatile T*  dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, T* dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
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
      //if (sz > 0) std::memset(dest,0,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(T* dest, std::size_t sz, const T& v) {
      //if (sz > 0) std::memset(dest,v,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(volatile T* dest, std::size_t sz, const T& v) {
      //if (sz > 0) std::memset(dest,v,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(std::size_t sz, const T& x = T(0.0)) {
      T* m = 0;
      if (sz > 0) {
        m = static_cast<T* >(operator new(sz*sizeof(T)));
        //std::memset(m,x,sz*sizeof(T));
        for (std::size_t i=0; i<sz; ++i)
          m[i] = x;
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
        m = static_cast<T* >(operator new(sz*sizeof(T)));
        for (std::size_t i=0; i<sz; ++i)
          m[i] = src[i];
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const volatile T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
        m = static_cast<T* >(operator new(sz*sizeof(T)));
        for (std::size_t i=0; i<sz; ++i)
          m[i] = src[i];
      }
      return m;
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(T* m, std::size_t sz) {
      if (sz > 0) operator delete((void*) m);
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(volatile T* m, std::size_t sz) {
      if (sz > 0) operator delete((void*) m);
    }
  };

  /*!
   * \brief Dynamic array allocation class that works for any type
   */
  template <typename T, typename device_t>
  struct DynArrayTraits<T, device_t, false> {

    typedef T value_type;
    typedef device_t execution_space;

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

    //! Get memory for new array of length \c sz and fill with zeros
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(std::size_t sz, const T& x = T(0.0)) {
      T* m = 0;
      if (sz > 0) {
        m = static_cast<T* >(operator new(sz*sizeof(T)));
        T* p = m;
        for (std::size_t i=0; i<sz; ++i)
          new (p++) T(x);
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
        m = static_cast<T* >(operator new(sz*sizeof(T)));
        T* p = m;
        for (std::size_t i=0; i<sz; ++i)
          new (p++) T(*(src++));
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const volatile T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
        m = static_cast<T* >(operator new(sz*sizeof(T)));
        T* p = m;
        for (std::size_t i=0; i<sz; ++i)
          new (p++) T(*(src++));
      }
      return m;
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(T* m, std::size_t sz) {
      T* e = m+sz;
      for (T* b = m; b!=e; b++)
        b->~T();
      operator delete((void*) m);
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(volatile T* m, std::size_t sz) {
      T* e = m+sz;
      for (T* b = m; b!=e; b++)
        b->~T();
      operator delete((void*) m);
    }
  };

#if defined(KOKKOS_ENABLE_CUDA)

  /*!
   * \brief Dynamic array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   *
   * Specialized for Cuda.
   */
  template <typename T>
  struct DynArrayTraits<T,Kokkos::Cuda,true> {

    typedef T value_type;
    typedef Kokkos::Cuda execution_space;

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, volatile T*  dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const volatile T* src, T*  dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, volatile T*  dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOS_INLINE_FUNCTION
    void copy(const T* src, T* dest, std::size_t sz) {
      //if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = *(src++);
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
      //if (sz > 0) std::memset(dest,0,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = T(0.);
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(T* dest, std::size_t sz, const T& v) {
      //if (sz > 0) std::memset(dest,v,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

    //! Fill array \c dest of length \c sz with value \c v
    static
    KOKKOS_INLINE_FUNCTION
    void fill(volatile T* dest, std::size_t sz, const T& v) {
      //if (sz > 0) std::memset(dest,v,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
        *(dest++) = v;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(std::size_t sz, const T& x = T(0.0)) {
      T* m = 0;
      if (sz > 0) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
        m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
        //std::memset(m,x,sz*sizeof(T));
        for (std::size_t i=0; i<sz; ++i)
          m[i] = x;
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
        m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
        for (std::size_t i=0; i<sz; ++i)
          m[i] = src[i];
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const volatile T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
        m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
        for (std::size_t i=0; i<sz; ++i)
          m[i] = src[i];
      }
      return m;
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(T* m, std::size_t sz) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
      cudaFree(m);
#else
      if (sz > 0) operator delete((void*) m);
#endif
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(volatile T* m, std::size_t sz) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
      cudaFree(m);
#else
      if (sz > 0) operator delete((void*) m);
#endif
    }
  };

  /*!
   * \brief Dynamic array allocation class that works for any type
   *
   * Specialized for Cuda
   */
  template <typename T>
  struct DynArrayTraits<T, Kokkos::Cuda, false> {

    typedef T value_type;
    typedef Kokkos::Cuda execution_space;

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

    //! Get memory for new array of length \c sz and fill with zeros
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(std::size_t sz, const T& x = T(0.0)) {
      T* m = 0;
      if (sz > 0) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
        m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
        T* p = m;
        for (std::size_t i=0; i<sz; ++i)
          new (p++) T(x);
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
        m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
        T* p = m;
        for (std::size_t i=0; i<sz; ++i)
          new (p++) T(*(src++));
      }
      return m;
    }

    /*!
     * \brief Get memory for new array of length \c sz and fill with
     * entries from \c src
     */
    static
    KOKKOS_INLINE_FUNCTION
    T* get_and_fill(const volatile T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_ENABLE_CUDA_UVM) && !defined( __CUDA_ARCH__ )
        cudaMallocManaged( (void**) &m, sz*sizeof(T), cudaMemAttachGlobal );
#else
        m = static_cast<T* >(operator new(sz*sizeof(T)));
#endif
        T* p = m;
        for (std::size_t i=0; i<sz; ++i)
          new (p++) T(*(src++));
      }
      return m;
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(T* m, std::size_t sz) {
      T* e = m+sz;
      for (T* b = m; b!=e; b++)
        b->~T();
      operator delete((void*) m);
    }

    //! Destroy array elements and release memory
    static
    KOKKOS_INLINE_FUNCTION
    void destroy_and_release(volatile T* m, std::size_t sz) {
      T* e = m+sz;
      for (T* b = m; b!=e; b++)
        b->~T();
      operator delete((void*) m);
    }
  };

#endif

} // namespace Stokhos

#endif // STOKHOS_DYN_ARRAY_TRAITS_HPP
