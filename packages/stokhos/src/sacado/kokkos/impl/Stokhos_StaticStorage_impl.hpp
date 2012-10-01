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

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Stokhos_StaticStorage_impl.hpp> without macros defined"

#else

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t, int Num>
  class StaticStorage<ordinal_t, value_t, Num, KOKKOSARRAY_MACRO_DEVICE> {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef KOKKOSARRAY_MACRO_DEVICE node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::StaticArrayTraits<value_type,node_type> ss;

    //! Turn StaticStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef StaticStorage<ord_t,val_t,Num,node_type> type;
    };

    //! Constructor
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    StaticStorage(const ordinal_type& sz,
		  const value_type& x = value_type(0.0)) : sz_(sz) { 
      ss::fill(coeff_, sz_, x); 
    }

    //! Copy constructor
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    StaticStorage(const StaticStorage& s) : sz_(s.sz_) {
      ss::copy(s.coeff_, coeff_, sz_);
    }

    //! Destructor
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    ~StaticStorage() {}

    //! Assignment operator
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    StaticStorage& operator=(const StaticStorage& s) {
      sz_ = s.sz_;
      ss::copy(s.coeff_, coeff_, sz_);
      return *this;
    }

    //! Initialize values to a constant value
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_reference v) { 
      ss::fill(coeff_, sz_, v); 
    }

    //! Initialize values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
      	ss::copy(v, coeff_, sz_);
      else
      	ss::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void load(pointer v) { 
      ss::copy(coeff_, v, sz_); 
    }

    //! Resize to new size (values are preserved)
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void resize(const ordinal_type& sz) { 
      if (sz > sz_)
	ss::fill(coeff_+sz_, sz-sz_, value_type(0.0));
      sz_ = sz; 
    }

    //! Reset storage to given array, size, and stride
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) {}

    //! Return size
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    ordinal_type size() const { return sz_; }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference operator[] (const ordinal_type& i) const { 
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference operator[] (const ordinal_type& i) { return coeff_[i]; }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference getCoeff() { return coeff_[i]; }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference getCoeff() const { return coeff_[i]; }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    pointer coeff() { return coeff_; }

  private:

    //! Coefficient values
    value_type coeff_[Num];

    //! Size of array used
    ordinal_type sz_;

  };

}

#endif // STOKHOS_STANDARD_STORAGE_HPP
