// $Id$
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_DYNAMIC_STRIDED_STORAGE_HPP
#define STOKHOS_DYNAMIC_STRIDED_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "KokkosArray_Macros.hpp"

namespace Stokhos {

  template <typename ordinal_t, typename value_t, typename node_t>
  class DynamicStridedStorage {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
    static const bool supports_reset = true;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef node_t node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::DynArrayTraits<value_type,node_type> ds;

    //! Turn DynamicStridedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef DynamicStridedStorage<ord_t,val_t,node_type> type;
    };

    //! Constructor
    KOKKOSARRAY_INLINE_FUNCTION
    DynamicStridedStorage(const ordinal_type& sz,
			  const value_type& x = value_type(0.0)) : 
      sz_(sz), stride_(1), is_owned_(true) {
      coeff_ = ds::get_and_fill(sz_, x);
    }

    //! Constructor
    KOKKOSARRAY_INLINE_FUNCTION
    DynamicStridedStorage(const DynamicStridedStorage& s) : 
    sz_(s.sz_), stride_(1), is_owned_(true) {
      if (s.stride_ == 1)
	coeff_ = ds::get_and_fill(s.coeff_, sz_);
      else {
	coeff_ = ds::get_and_fill(sz_);
	for (ordinal_type i=0; i<sz_; ++i)
	  coeff_[i] = s[i];
      }
    }

    //! Destructor
    KOKKOSARRAY_INLINE_FUNCTION
    ~DynamicStridedStorage() {
      if (is_owned_) ds::destroy_and_release(coeff_, sz_*stride_);
    }

    //! Assignment operator
    KOKKOSARRAY_INLINE_FUNCTION
    DynamicStridedStorage& operator=(const DynamicStridedStorage& s) {
      if (&s != this) { 
	if (s.sz_ != sz_) {
	  if (is_owned_)
	    ds::destroy_and_release(coeff_, sz_*stride_);
	  if (s.stride_ == 1)
	    coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
	  else {
	    coeff_ = ds::get_and_fill(s.sz_);
	    for (ordinal_type i=0; i<s.sz_; ++i)
	      coeff_[i] = s[i];
	  }
	  sz_ = s.sz_;
	  stride_ = 1;
	  is_owned_ = true;
	}
	else {
	  if (stride_ == 1 and s.stride_ == 1)
	    ds::copy(s.coeff_, coeff_, sz_);
	  else
	    for (ordinal_type i=0; i<s.sz_; ++i)
	      coeff_[i*stride_] = s[i];
	}
      }
      return *this;
    }

    //! Initialize values to a constant value
    KOKKOSARRAY_INLINE_FUNCTION
    void init(const_reference v) { 
      if (stride_ == 1)
	ds::fill(coeff_, sz_, v); 
      else
	for (ordinal_type i=0; i<sz_; ++i)
	  coeff_[i*stride_] = v;
    }

    //! Initialize values to an array of values
    KOKKOSARRAY_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      ordinal_type my_sz = sz;
      if (sz == 0)
	my_sz = sz_;
      if (stride_ == 1)
	ds::copy(v, coeff_, my_sz);
      else
	for (ordinal_type i=0; i<my_sz; ++i)
	  coeff_[i*stride_] = v[i];
    }

    //! Load values to an array of values
    KOKKOSARRAY_INLINE_FUNCTION
    void load(pointer v) {
      if (stride_ == 1)
	ds::copy(coeff_, v, sz_); 
      for (ordinal_type i=0; i<sz_; ++i)
	coeff_[i*stride_] = v[i];
    }

    //! Resize to new size (values are preserved)
    KOKKOSARRAY_INLINE_FUNCTION
    void resize(const ordinal_type& sz) { 
      if (sz != sz_) {
	value_type *coeff_new = ds::get_and_fill(sz);
	ordinal_type my_sz = sz_;
	if (sz_ > sz)
	  my_sz = sz;
	if (stride_ == 1)
	  ds::copy(coeff_, coeff_new, my_sz);
	else
	  for (ordinal_type i=0; i<my_sz; ++i)
	    coeff_new[i] = coeff_[i*stride_];
	if (is_owned_)
	  ds::destroy_and_release(coeff_, sz_*stride_);
	coeff_ = coeff_new;
	sz_ = sz;
	stride_ = 1;
	is_owned_ = true;
      }
    }

    //! Reset storage to given array, size, and stride
    KOKKOSARRAY_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) { 
      if (is_owned_)
	ds::destroy_and_release(coeff_, sz_*stride_);
      coeff_ = v;
      sz_ = sz;
      stride_ = stride;
      is_owned_ = owned;
    }

    //! Return size
    KOKKOSARRAY_INLINE_FUNCTION
    ordinal_type size() const { return sz_; }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_INLINE_FUNCTION
    const_reference operator[] (const ordinal_type& i) const {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_INLINE_FUNCTION
    reference operator[] (const ordinal_type& i) {
      return coeff_[i*stride_];
    }

    template <int i>
    KOKKOSARRAY_INLINE_FUNCTION
    reference getCoeff() { return coeff_[i]; }

    template <int i>
    KOKKOSARRAY_INLINE_FUNCTION
    const_reference getCoeff() const { return coeff_[i]; }

    //! Get coefficients
    KOKKOSARRAY_INLINE_FUNCTION
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    KOKKOSARRAY_INLINE_FUNCTION
    pointer coeff() { return coeff_; }

  private:

    //! Coefficient values
    pointer coeff_;

    //! Size of array used
    ordinal_type sz_;

    //! Stride of array
    ordinal_type stride_;

    //! Do we own the array
    bool is_owned_;

  };

}

#endif // STOKHOS_DYNAMIC_STORAGE_HPP
