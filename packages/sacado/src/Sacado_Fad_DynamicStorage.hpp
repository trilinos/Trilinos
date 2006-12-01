// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef SACADO_FAD_DYNAMICSTORAGE_HPP
#define SACADO_FAD_DYNAMICSTORAGE_HPP

#include <new>

#include "Sacado_Traits.hpp"

namespace Sacado {

  namespace Fad {

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
	  new (p++) T(0.);
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
	memset(m,0,sz*sizeof(T));
	return m;
      }

      /*! 
       * \brief Get memory for new array of length \c sz and fill with 
       * entries from \c src
       */
      static inline T* get_and_fill(const T* src, int sz) {
	T* m = new T[sz];
	for (int i=0; i<sz; ++i)
	  m[i] = src[i];
	return m;
      }

      //! Copy array from \c src to \c dest of length \c sz
      static inline void copy(const T* src, T* dest, int sz) {
	memcpy(dest,src,sz*sizeof(T));
      }

      //! Zero out array \c dest of length \c sz
      static inline void zero(T* dest, int sz) {
	memset(dest,0,sz*sizeof(T));
      }
    };

    //! Derivative array storage class using dynamic memory allocation
    template <typename T> 
    class DynamicStorage {

    public:

      //! Default constructor
      DynamicStorage() : sz_(0), len_(0), dx_(NULL) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      DynamicStorage(const int sz) : sz_(sz), len_(sz) {
	dx_ = ds_array<T>::get_and_fill(sz_);
      }

      //! Copy constructor
      DynamicStorage(const DynamicStorage& x) : sz_(x.sz_), len_(x.sz_) {
	dx_ = ds_array<T>::get_and_fill(x.dx_, sz_);
      }
      
      //! Destructor
      ~DynamicStorage() {
	if (len_ != 0)
	  operator delete((void*)dx_);
      }

      //! Assignment
      DynamicStorage& operator=(const DynamicStorage& x) { 
	if (sz_ != x.sz_) {
	  sz_ = x.sz_;
	  if (x.sz_ > len_) {
	    if (len_ != 0)
	      operator delete((void*)dx_);
	    len_ = x.sz_;
	    dx_ = ds_array<T>::get_and_fill(x.dx_, sz_);
	  }
	  else 
	    ds_array<T>::copy(x.dx_, dx_, sz_);
	}
	else 
	  ds_array<T>::copy(x.dx_, dx_, sz_);

	return *this; 
      } 

      //! Returns number of derivative components
      int size() const { return sz_;}

      //! Resize the derivative array to sz
      void resize(int sz) { 
	if (sz > len_) {
	  if (len_ != 0)
	    operator delete((void*)dx_);
	  dx_ = static_cast<T* >(operator new(sz*sizeof(T)));
	  len_ = sz;
	}
	sz_ = sz;
      }

      //! Zero out derivative array
      void zero() { 
	ds_array<T>::zero(dx_, sz_);
      }

    public:

      //! Derivative array size
      int sz_;

      //! Derivative array length
      int len_;

      //! Derivative array
      T* dx_;

    }; // class DynamicStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_DYNAMICSTORAGE_HPP
