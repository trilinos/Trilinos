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

#ifndef STOKHOS_DYNAMIC_STORAGE_HPP
#define STOKHOS_DYNAMIC_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "KokkosArray_Macros.hpp"

namespace Stokhos {

  template <typename ordinal_t, typename value_t, typename node_t>
  class DynamicStorage {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef node_t node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::DynArrayTraits<value_type,node_type> ds;

    //! Turn DynamicStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef DynamicStorage<ord_t,val_t,node_type> type;
    };

    //! Constructor
    KOKKOSARRAY_INLINE_FUNCTION
    DynamicStorage(const ordinal_type& sz,
		   const value_type& x = value_type(0.0)) : sz_(sz) {
      coeff_ = ds::get_and_fill(sz_, x);
    }

    //! Constructor
    KOKKOSARRAY_INLINE_FUNCTION
    DynamicStorage(const DynamicStorage& s) : sz_(s.sz_) {
      coeff_ = ds::get_and_fill(s.coeff_, sz_);
    }

    //! Destructor
    KOKKOSARRAY_INLINE_FUNCTION
    ~DynamicStorage() {
      ds::destroy_and_release(coeff_, sz_);
    }

    //! Assignment operator
    KOKKOSARRAY_INLINE_FUNCTION
    DynamicStorage& operator=(const DynamicStorage& s) {
      if (&s != this) { 
	if (s.sz_ != sz_) {
	  ds::destroy_and_release(coeff_, sz_);
	  coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
	  sz_ = s.sz_;
	}
	else
	  ds::copy(s.coeff_, coeff_, sz_);
      }
      return *this;
    }

    //! Initialize values to a constant value
    KOKKOSARRAY_INLINE_FUNCTION
    void init(const_reference v) { 
      ds::fill(coeff_, sz_, v); 
    }

    //! Initialize values to an array of values
    KOKKOSARRAY_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
      	ds::copy(v, coeff_, sz_);
      else
      	ds::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOSARRAY_INLINE_FUNCTION
    void load(pointer v) {
      ds::copy(coeff_, v, sz_); 
    }

    //! Resize to new size (values are preserved)
    KOKKOSARRAY_INLINE_FUNCTION
    void resize(const ordinal_type& sz) { 
      if (sz != sz_) {
	value_type *coeff_new = ds::get_and_fill(sz);
	if (sz > sz_)
	  ds::copy(coeff_, coeff_new, sz_);
	else
	  ds::copy(coeff_, coeff_new, sz);
	ds::destroy_and_release(coeff_, sz_);
	coeff_ = coeff_new;
	sz_ = sz;
      }
    }

    //! Reset storage to given array, size, and stride
    KOKKOSARRAY_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) {}

    //! Return size
    KOKKOSARRAY_INLINE_FUNCTION
    ordinal_type size() const { return sz_; }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_INLINE_FUNCTION
    const_reference operator[] (const ordinal_type& i) const {
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_INLINE_FUNCTION
    reference operator[] (const ordinal_type& i) {
      return coeff_[i];
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

  };

}

#endif // STOKHOS_DYNAMIC_STORAGE_HPP
