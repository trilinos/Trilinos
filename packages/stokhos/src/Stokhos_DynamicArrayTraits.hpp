// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_DYNAMICARRAYTRAITS_HPP
#define STOKHOS_DYNAMICARRAYTRAITS_HPP

#include <new>
#include <cstring>

namespace Stokhos {

  //! Base template specification for %IsScalarType
  /*!
   * The %IsScalarType classes provide a mechanism for computing the 
   * determining whether a type is a scalar type (float, double, etc...)
   */
  template <typename T> struct IsScalarType {
    static const bool value = false;
  };

  //! Specialization of above classes to built-in types
#define STOKHOS_BUILTIN_SPECIALIZATION(t)                  \
  template <> struct IsScalarType< t > {	          \
    static const bool value = true;	       		  \
  };

  STOKHOS_BUILTIN_SPECIALIZATION(float)
  STOKHOS_BUILTIN_SPECIALIZATION(double)
  STOKHOS_BUILTIN_SPECIALIZATION(int)
  STOKHOS_BUILTIN_SPECIALIZATION(long)

#undef STOKHOS_BUILTIN_SPECIALIZATION
  

  /*!
   * \brief Dynamic array allocation class that works for any type
   */
  template <typename T, bool isScalar = IsScalarType<T>::value>
  struct ds_array {

    //! Get memory for new array of length \c sz and fill with zeros
    static inline T* get_and_fill(int sz) {
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
      T* p = m;
      for (int i=0; i<sz; ++i)
	new (p++) T(0.0);
      return m;
    }

    /*! 
     * \brief Get memory for new array of length \c sz and fill with 
     * entries from \c src
     */
    static inline T* get_and_fill(const T* src, int sz) {
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
      T* p = m; 
      for (int i=0; i<sz; ++i)
	new (p++) T(*(src++));
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T*  dest, int sz) {
      for (int i=0; i<sz; ++i)
	*(dest++) = *(src++);
    }

    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, int sz) {
      for (int i=0; i<sz; ++i)
	*(dest++) = T(0.);
    }

    //! Destroy array elements and release memory
    static inline void destroy_and_release(T* m, int sz) {
      T* e = m+sz;
      for (T* b = m; b!=e; b++)
	b->~T();
      operator delete((void*) m);
    }
  };

  /*!
   * \brief Dynamic array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T>
  struct ds_array<T,true> {

    //! Get memory for new array of length \c sz and fill with zeros
    static inline T* get_and_fill(int sz) {
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
      std::memset(m,0,sz*sizeof(T));
      return m;
    }

    /*! 
     * \brief Get memory for new array of length \c sz and fill with 
     * entries from \c src
     */
    static inline T* get_and_fill(const T* src, int sz) {
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
      for (int i=0; i<sz; ++i)
	m[i] = src[i];
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T* dest, int sz) {
      std::memcpy(dest,src,sz*sizeof(T));
    }

    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, int sz) {
      std::memset(dest,0,sz*sizeof(T));
    }

    //! Destroy array elements and release memory
    static inline void destroy_and_release(T* m, int sz) {
      operator delete((void*) m);
      }
  };

} // namespace Stokhos

#endif // STOKHOS_DYNAMICARRAYTRAITS_HPP
