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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
    StaticStandardStorage(const ordinal_type& sz_,
			  const value_type& x = value_type(0.0)) : sz(sz_) { 
      std::fill(coeff_, coeff_+sz, x); 
    }

    //! Copy constructor
    StaticStandardStorage(const StaticStandardStorage& s) : sz(s.sz) {
      ss::copy(s.coeff_, coeff_, sz);
    }

    //! Destructor
    ~StaticStandardStorage() {}

    //! Assignment operator
    StaticStandardStorage& operator=(const StaticStandardStorage& s) {
      sz = s.sz;
      ss::copy(s.coeff_, coeff_, sz);
      return *this;
    }

    //! Initialize values to a constant value
    void init(const_reference v) { 
      std::fill(coeff_, coeff_+sz, v); 
    }

    //! Initialize values to an array of values
    void init(const_pointer v, const ordinal_type& sz_ = 0) {
      if (sz_ == 0)
      	ss::copy(v, coeff_, sz);
      else
      	ss::copy(v, coeff_, sz_);
    }

    //! Load values to an array of values
    void load(pointer v) { 
      ss::copy(coeff_, v, sz); 
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
