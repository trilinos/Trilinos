// $Id$ 
// $Source$ 
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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

namespace Sacado {

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

    /*! 
     * \brief Get memory for new array of length \c sz and fill with 
     * entries from \c src
     */
    static inline T* strided_get_and_fill(const T* src, int stride, int sz) {
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
      T* p = m; 
      for (int i=0; i<sz; ++i) {
	new (p++) T(*(src));
	src += stride;
      }
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T*  dest, int sz) {
      for (int i=0; i<sz; ++i)
	*(dest++) = *(src++);
    }

    //! Copy array from \c src to \c dest of length \c sz
    static inline void strided_copy(const T* src, int src_stride, 
				    T* dest, int dest_stride, int sz) {
      for (int i=0; i<sz; ++i) {
	*(dest) = *(src);
	dest += dest_stride;
	src += src_stride;
      }
    }

    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, int sz) {
      for (int i=0; i<sz; ++i)
	*(dest++) = T(0.);
    }

    //! Zero out array \c dest of length \c sz
    static inline void strided_zero(T* dest, int stride, int sz) {
      for (int i=0; i<sz; ++i) {
	*(dest) = T(0.);
	dest += stride;
      }
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
      if (sz > 0) {
	T* m = static_cast<T* >(operator new(sz*sizeof(T)));
	std::memset(m,0,sz*sizeof(T));
	return m;
      }
      else
	return NULL;
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

    /*! 
     * \brief Get memory for new array of length \c sz and fill with 
     * entries from \c src
     */
    static inline T* strided_get_and_fill(const T* src, int stride, int sz) {
      T* m = static_cast<T* >(operator new(sz*sizeof(T)));
      for (int i=0; i<sz; ++i)
	m[i] = src[i*stride];
      return m;
    }

    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T* dest, int sz) {
      std::memcpy(dest,src,sz*sizeof(T));
    }

    //! Copy array from \c src to \c dest of length \c sz
    static inline void strided_copy(const T* src, int src_stride, 
				    T* dest, int dest_stride, int sz) {
      for (int i=0; i<sz; ++i) {
	*(dest) = *(src);
	dest += dest_stride;
	src += src_stride;
      }
    }

    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, int sz) {
      if (sz > 0)
	std::memset(dest,0,sz*sizeof(T));
    }

    //! Zero out array \c dest of length \c sz
    static inline void strided_zero(T* dest, int stride, int sz) {
      for (int i=0; i<sz; ++i) {
	*(dest) = T(0.);
	dest += stride;
      }
    }

    //! Destroy array elements and release memory
    static inline void destroy_and_release(T* m, int sz) {
      operator delete((void*) m);
    }
  };

} // namespace Sacado

#endif // SACADO_DYNAMICARRAY_HPP
