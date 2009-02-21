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

#ifndef SACADO_FAD_VECTORDYNAMICSTORAGE_HPP
#define SACADO_FAD_VECTORDYNAMICSTORAGE_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_DynamicArrayTraits.hpp"

namespace Sacado {

  namespace Fad {

    //! Derivative array storage class using dynamic memory allocation
    template <typename T, typename S = T> 
    class VectorDynamicStorage {

    public:

      //! Default constructor
      VectorDynamicStorage(const T & x) : 
	owns_mem(true), sz_(0), len_(0), val_(new T(x)), dx_(NULL) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      VectorDynamicStorage(const int sz, const T & x) : 
	owns_mem(true), sz_(sz), len_(sz), val_(new T(x)) {
	dx_ = ds_array<S>::get_and_fill(sz_);
      }

      //! Constructor with supplied memory
      VectorDynamicStorage(const int sz, T* x, S* dx, bool zero_out) : 
	owns_mem(false), sz_(sz), len_(sz), val_(x), dx_(dx) {
	if (zero_out)
	  zero(dx_, sz_);
      }

      //! Copy constructor
      VectorDynamicStorage(const VectorDynamicStorage& x) : 
	owns_mem(true), sz_(x.sz_), len_(x.sz_), val_(new T(*x.val_))  {
	dx_ = ds_array<S>::get_and_fill(x.dx_, sz_);
      }
      
      //! Destructor
      ~VectorDynamicStorage() {
	if (owns_mem) {
	  delete val_;
	  if (len_ != 0)
	    ds_array<S>::destroy_and_release(dx_, len_);
	}
      }

      //! Assignment
      VectorDynamicStorage& operator=(const VectorDynamicStorage& x) { 
	*val_ = *x.val_;
	if (sz_ != x.sz_) {
	  sz_ = x.sz_;
	  if (x.sz_ > len_) {
	    if (!owns_mem)
	      throw "Can\'t resize beyond original size when memory isn't owned!";
	    if (len_ != 0)
	      ds_array<S>::destroy_and_release(dx_, len_);
	    len_ = x.sz_;
	    dx_ = ds_array<S>::get_and_fill(x.dx_, sz_);
	  }
	  else 
	    ds_array<S>::copy(x.dx_, dx_, sz_);
	}
	else 
	  ds_array<S>::copy(x.dx_, dx_, sz_);

	return *this; 
      } 

      //! Returns number of derivative components
      int size() const { return sz_;}

      //! Resize the derivative array to sz
      void resize(int sz) { 
	if (sz > len_) {
	  if (!owns_mem)
	      throw "Can\'t resize beyond original size when memory isn't owned!";
	  if (len_ != 0)
	    ds_array<S>::destroy_and_release(dx_, len_);
	  dx_ = ds_array<S>::get_and_fill(sz);
	  len_ = sz;
	}
	sz_ = sz;
      }

      //! Zero out derivative array
      void zero() { 
	ds_array<S>::zero(dx_, sz_);
      }

      //! Set value/derivative array memory
      void setMemory(int sz, T* x, S* dx) {

	// Destroy old memory
	if (owns_mem) {
	  delete val_;
	  if (len_ != 0)
	    ds_array<S>::destroy_and_release(dx_, len_);
	}

	// Set new values
	owns_mem = false;
	sz_ = sz;
	len_ = sz;
	val_ = x;
	dx_ = dx;
      }

    public:

      //! Do we own the val/dx storage
      bool owns_mem;

      //! Derivative array size
      int sz_;

      //! Derivative array length
      int len_;

      //! Value
      T* val_;

      //! Derivative array
      S* dx_;

    }; // class VectorDynamicStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_VECTORDYNAMICSTORAGE_HPP
