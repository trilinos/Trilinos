// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_STATIC_FIXED_STANDARD_STORAGE_HPP
#define STOKHOS_STATIC_FIXED_STANDARD_STORAGE_HPP

#include "Sacado_StaticArrayTraits.hpp"
#include <algorithm>

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_type, typename value_type, int Num>
  class StaticFixedStandardStorage {
  public:

    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Sacado::ss_array<value_type> ss;

    //! Turn StaticFixedStandardStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef StaticFixedStandardStorage<ord_t,val_t,Num> type;
    };

    //! Constructor
    StaticFixedStandardStorage(const ordinal_type& sz,
			       const value_type& x = value_type(0.0)) { 
      std::fill(coeff_, coeff_+Num, x); 
    }

    //! Copy constructor
    StaticFixedStandardStorage(const StaticFixedStandardStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Destructor
    ~StaticFixedStandardStorage() {}

    //! Assignment operator
    StaticFixedStandardStorage& operator=(const StaticFixedStandardStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
      return *this;
    }

    //! Initialize values to a constant value
    void init(const_reference v) { 
      std::fill(coeff_, coeff_+Num, v); 
    }

    //! Initialize values to an array of values
    void init(const_pointer v, const ordinal_type& sz_ = 0) {
      if (sz_ == 0)
      	ss::copy(v, coeff_, Num);
      else
      	ss::copy(v, coeff_, sz_);
    }

    //! Load values to an array of values
    void load(pointer v) { 
      ss::copy(coeff_, v, Num); 
    }

    //! Resize to new size (values are preserved)
    void resize(const ordinal_type& sz) {}

    //! Return size
    static ordinal_type size() { return Num; }

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

  };

}

#endif // STOKHOS_STANDARD_STORAGE_HPP
