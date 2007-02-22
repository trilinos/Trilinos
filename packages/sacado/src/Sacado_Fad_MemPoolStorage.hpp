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

#ifndef SACADO_FAD_MEMPOOLSTORAGE_HPP
#define SACADO_FAD_MEMPOOLSTORAGE_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_Traits.hpp"
#include "Sacado_Fad_MemPool.hpp"

namespace Sacado {

  namespace Fad {

    /*!
     * \brief Dynamic array allocation class that works for any type
     */
    template <typename T, bool isScalar = IsScalarType<T>::value>
    struct mp_array {

      //! Get memory for new array of length \c sz and fill with zeros
      static inline T* get_and_fill(int sz, MemPool* pool) {
	if (sz) {
	  T* m = static_cast<T*>(pool->alloc());
	  T* p = m;
	  for (int i=0; i<sz; ++i)
	    new (p++) T(0.);
	  return m;
	}
	else
	  return NULL;
      }

      /*! 
       * \brief Get memory for new array of length \c sz and fill with 
       * entries from \c src
       */
      static inline T* get_and_fill(const T* src, int sz, MemPool* pool) {
	if (sz) {
	  T* m = static_cast<T*>(pool->alloc());
	  T* p = m; 
	  for (int i=0; i<sz; ++i)
	    new (p++) T(*(src++));
	  return m;
	}
	else
	  return NULL;
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
      static inline void destroy_and_release(T* m, int sz, MemPool* pool) {
	T* e = m+sz;
	for (T* b = m; b!=e; b++)
	  b->~T();
	pool->free((void*) m);
      }
    };

    /*!
     * \brief Dynamic array allocation class that is specialized for scalar
     * i.e., fundamental or built-in types (float, double, etc...).
     */
    template <typename T>
    struct mp_array<T,true> {

      //! Get memory for new array of length \c sz and fill with zeros
      static inline T* get_and_fill(int sz, MemPool* pool) {
	if (sz) {
	  T* m = static_cast<T*>(pool->alloc());
	  memset(m,0,sz*sizeof(T));
	  return m;
	}
	else
	  return NULL;
      }

      /*! 
       * \brief Get memory for new array of length \c sz and fill with 
       * entries from \c src
       */
      static inline T* get_and_fill(const T* src, int sz, MemPool* pool) {
	if (sz) {
	  T* m = static_cast<T*>(pool->alloc());
	  T* p = m; 
	  for (int i=0; i<sz; ++i)
	    new (p++) T(*(src++));
	  return m;
	}
	else
	  return NULL;
      }

      //! Copy array from \c src to \c dest of length \c sz
      static inline void copy(const T* src, T* dest, int sz) {
	memcpy(dest,src,sz*sizeof(T));
      }

      //! Zero out array \c dest of length \c sz
      static inline void zero(T* dest, int sz) {
	memset(dest,0,sz*sizeof(T));
      }

      //! Destroy array elements and release memory
      static inline void destroy_and_release(T* m, int sz, MemPool* pool) {
	pool->free((void*) m);
      }
    };

    //! Derivative array storage class using dynamic memory allocation
    template <typename T> 
    class MemPoolStorage {

    public:

      //! Default constructor
      MemPoolStorage(const T & x) : 
	val_(x), sz_(0), len_(0), dx_(NULL), myPool_(defaultPool_) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      MemPoolStorage(const int sz, const T & x) : 
	val_(x), sz_(sz), len_(sz), myPool_(defaultPool_) {
	dx_ = mp_array<T>::get_and_fill(sz_, myPool_);
      }

      //! Copy constructor
      MemPoolStorage(const MemPoolStorage& x) : 
	val_(x.val_), sz_(x.sz_), len_(x.sz_), myPool_(x.myPool_) {
	dx_ = mp_array<T>::get_and_fill(x.dx_, sz_, myPool_);
      }
      
      //! Destructor
      ~MemPoolStorage() {
	if (len_ != 0)
	  mp_array<T>::destroy_and_release(dx_, len_, myPool_);
      }

      //! Assignment
      MemPoolStorage& operator=(const MemPoolStorage& x) { 
	val_ = x.val_;
	if (sz_ != x.sz_) {
	  sz_ = x.sz_;
	  if (x.sz_ > len_) {
	    if (len_ != 0)
	      mp_array<T>::destroy_and_release(dx_, len_, myPool_);
	    len_ = x.sz_;
	    myPool_ = x.myPool_;
	    dx_ = mp_array<T>::get_and_fill(x.dx_, sz_, myPool_);
	  }
	  else 
	    mp_array<T>::copy(x.dx_, dx_, sz_);
	}
	else 
	  mp_array<T>::copy(x.dx_, dx_, sz_);

	return *this; 
      } 

      //! Returns number of derivative components
      int size() const { return sz_;}

      //! Resize the derivative array to sz
      void resize(int sz) { 
	if (sz > len_) {
	   if (len_ != 0)
	    mp_array<T>::destroy_and_release(dx_, len_, myPool_);
	  myPool_ = defaultPool_;
	  dx_ = mp_array<T>::get_and_fill(sz, myPool_);
	  len_ = sz;
	}
	sz_ = sz;
      }

      //! Zero out derivative array
      void zero() { 
	mp_array<T>::zero(dx_, sz_);
      }

    public:

      //! Value
      T val_;

      //! Derivative array size
      int sz_;

      //! Derivative array length
      int len_;

      //! Derivative array
      T* dx_;

      //! Default memory pool
      static MemPool* defaultPool_;

    protected:

      //! Memory pool
      MemPool* myPool_;

    }; // class MemPoolStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_MEMPOOLSTORAGE_HPP
