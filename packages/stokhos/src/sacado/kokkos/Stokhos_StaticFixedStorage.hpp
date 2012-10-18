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

#ifndef STOKHOS_STATIC_FIXED_STORAGE_HPP
#define STOKHOS_STATIC_FIXED_STORAGE_HPP

#include "Stokhos_StaticArrayTraits.hpp"

#include "KokkosArray_Macros.hpp"

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t, int Num, typename node_t>
  class StaticFixedStorage {
  public:

    static const bool is_static = true;
    static const int static_size = Num;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef node_t node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::StaticArrayTraits<value_type,node_type> ss;

    //! Turn StaticFixedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef StaticFixedStorage<ord_t,val_t,Num,node_type> type;
    };

    //! Constructor
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz,
		       const value_type& x = value_type(0.0)) { 
      ss::fill(coeff_, Num, x); 
    }

    //! Copy constructor
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    StaticFixedStorage(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Destructor
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    ~StaticFixedStorage() {}

    //! Assignment operator
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    StaticFixedStorage& operator=(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
      return *this;
    }

    //! Initialize values to a constant value
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    void init(const_reference v) { 
      ss::fill(coeff_, Num, v); 
    }

    //! Initialize values to an array of values
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
      	ss::copy(v, coeff_, Num);
      else
      	ss::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    void load(pointer v) { 
      ss::copy(coeff_, v, Num); 
    }

    //! Resize to new size (values are preserved)
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    void resize(const ordinal_type& sz) {}

    //! Reset storage to given array, size, and stride
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) {}

    //! Return size
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    static ordinal_type size() { return Num; }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    const_reference operator[] (const ordinal_type& i) const { 
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    reference operator[] (const ordinal_type& i) { return coeff_[i]; }

    template <int i>
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    reference getCoeff() { return coeff_[i]; }

    template <int i>
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    const_reference getCoeff() const { return coeff_[i]; }

    //! Get coefficients
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    pointer coeff() { return coeff_; }

  private:

    //! Coefficient values
    value_type coeff_[Num];

  };

}

#endif // STOKHOS_STATIC_FIXED_STORAGE_HPP
