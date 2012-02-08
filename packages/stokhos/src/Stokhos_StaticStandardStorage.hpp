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

#ifndef STOKHOS_STATIC_STANDARD_STORAGE_HPP
#define STOKHOS_STATIC_STANDARD_STORAGE_HPP

#include "Sacado_StaticArrayTraits.hpp"
#include <algorithm>

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_type, typename value_type, int Num>
  class StaticStandardStorage {
  public:

    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Sacado::ss_array<value_type> ss;

    //! Turn StaticStandardStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef StaticStandardStorage<ord_t,val_t,Num> type;
    };

    //! Constructor
    StaticStandardStorage(const ordinal_type& sz_) : sz(sz_) { 
      //ss::zero(coeff_, sz); 
      std::fill(coeff_, coeff_+sz, value_type(0.0)); 
    }

    //! Copy constructor
    StaticStandardStorage(const StaticStandardStorage& s) : sz(s.sz) {
      //ss::copy(s.coeff_, coeff_, sz);
      std::copy(s.coeff_, s.coeff_+sz, coeff_);
    }

    //! Destructor
    ~StaticStandardStorage() {}

    //! Assignment operator
    StaticStandardStorage& operator=(const StaticStandardStorage& s) {
      sz = s.sz;
      //ss::copy(s.coeff_, coeff_, sz);
      std::copy(s.coeff_, s.coeff_+sz, coeff_);
      return *this;
    }

    //! Initialize values to a constant value
    void init(const_reference v) { std::fill(coeff_, coeff_+sz, v); }

    //! Initialize values to an array of values
    void init(const_pointer v, const ordinal_type& sz_ = 0) {
      // if (sz_ == 0)
      // 	ss::copy(v, coeff_, sz);
      // else
      // 	ss::copy(v, coeff_, sz_);

      if (sz_ == 0)
	std::copy(v, v+sz, coeff_);
      else
	std::copy(v, v+sz_, coeff_);
    }

    //! Load values to an array of values
    void load(pointer v) { 
      //ss::copy(coeff_, v, sz); 
      std::copy(coeff_, coeff_+sz, v);
    }

    //! Resize to new size (values are preserved)
    void resize(const ordinal_type& sz_) { 
      if (sz_ > sz)
	std::fill(coeff_+sz, coeff_+sz_, value_type(0.0));
      sz = sz_; 
    }

    //! Return size
    ordinal_type size() const { return sz; }

    //! Coefficient access (avoid if possible)
    const_reference operator[] (const ordinal_type& i) const { 
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    reference operator[] (const ordinal_type& i) { return coeff_[i]; }

    //! Get coefficients
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    pointer coeff() { return coeff_; }

  private:

    //! Coefficient values
    value_type coeff_[Num];

    //! Size of array used
    ordinal_type sz;

  };

}

#endif // STOKHOS_STANDARD_STORAGE_HPP
