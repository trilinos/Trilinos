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

#ifndef STOKHOS_STANDARD_STORAGE_HPP
#define STOKHOS_STANDARD_STORAGE_HPP

#include "Teuchos_Array.hpp"
#include <algorithm>

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class StandardStorage {
  public:

    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;

    //! Turn StandardStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef StandardStorage<ord_t,val_t> type;
    };

    //! Constructor
    StandardStorage(const ordinal_type& sz) : 
      coeff_(sz, value_type(0.0)) {}

    //! Destructor
    ~StandardStorage() {}

    //! Initialize values to a constant value
    void init(const_reference v) { 
      std::fill(coeff_.begin(), coeff_.end(), v); 
    }

    //! Initialize values to an array of values
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
	std::copy(v, v+coeff_.size(), coeff_.begin());
      else
	std::copy(v, v+sz, coeff_.begin());
    }

    //! Load values to an array of values
    void load(pointer v) {
      std::copy(coeff_.begin(), coeff_.end(), v);
    }

    //! Resize to new size (values are preserved)
    void resize(const ordinal_type& sz) { coeff_.resize(sz); }

    //! Return size
    ordinal_type size() const { return coeff_.size(); }

    //! Coefficient access (avoid if possible)
    const_reference operator[] (const ordinal_type& i) const {
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    reference operator[] (const ordinal_type& i) {
      return coeff_[i];
    }

    //! Get coefficients
    const_pointer coeff() const { return coeff_.getRawPtr(); }

    //! Get coefficients
    pointer coeff() { return coeff_.getRawPtr(); }

  private:

    //! Coefficient values
    Teuchos::Array<value_type> coeff_;

  };

}

#endif // STOKHOS_STANDARD_STORAGE_HPP
