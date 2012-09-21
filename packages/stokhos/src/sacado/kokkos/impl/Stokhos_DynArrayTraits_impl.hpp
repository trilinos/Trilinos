// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Stokhos_DynArrayTraits_impl.hpp> without macros defined"

#else

#include <new>
#include <cstring>

namespace Stokhos {

  /*!
   * \brief Dynamic array allocation class that works for any type
   */
  template <typename T>
  struct DynArrayTraits<T, KOKKOSARRAY_MACRO_DEVICE, false> {

    typedef T value_type;
    typedef KOKKOSARRAY_MACRO_DEVICE node_type;

    //! Fill array \c dest of length \c sz with value \c v
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void fill(T* dest, std::size_t sz, const T& v) {
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = v;
    }

    //! Copy array from \c src to \c dest of length \c sz
    static inline 
     KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void copy(const T* src, T*  dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = *(src++);
    }

    //! Zero out array \c dest of length \c sz
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void zero(T* dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = T(0.);
    }

    //! Get memory for new array of length \c sz and fill with zeros
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
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
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
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

    //! Destroy array elements and release memory
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void destroy_and_release(T* m, std::size_t sz) {
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
  struct DynArrayTraits<T,KOKKOSARRAY_MACRO_DEVICE,true> {

    typedef T value_type;
    typedef KOKKOSARRAY_MACRO_DEVICE node_type;

    //! Copy array from \c src to \c dest of length \c sz
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void copy(const T* src, T* dest, std::size_t sz) {
      if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
    }

    //! Zero out array \c dest of length \c sz
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void zero(T* dest, std::size_t sz) {
      if (sz > 0) std::memset(dest,0,sz*sizeof(T));
    }

    //! Fill array \c dest of length \c sz with value \c v
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void fill(T* dest, std::size_t sz, const T& v) {
      //if (sz > 0) std::memset(dest,v,sz*sizeof(T));
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = v;
    }

    //! Get memory for new array of length \c sz and fill with zeros
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
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
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    T* get_and_fill(const T* src, std::size_t sz) {
      T* m = 0;
      if (sz > 0) {
	m = static_cast<T* >(operator new(sz*sizeof(T)));
	for (std::size_t i=0; i<sz; ++i)
	  m[i] = src[i];
      }
      return m;
    }

    //! Destroy array elements and release memory
    static inline 
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void destroy_and_release(T* m, std::size_t sz) {
      if (sz > 0) operator delete((void*) m);
    }
  };

} // namespace Stokhos

#endif // STOKHOS_DYNAMICARRAYTRAITS_HPP
