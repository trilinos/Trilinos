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

#error "Including <Stokhos_LocalStorage_impl.hpp> without macros defined"

#else

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t>
  class LocalStorage<ordinal_t, value_t, 2, KOKKOSARRAY_MACRO_DEVICE> {
  public:

    static const int Num = 2;
    
    static const bool is_static = true;
    static const int static_size = Num;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef KOKKOSARRAY_MACRO_DEVICE node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;

    //! Turn LocalStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef LocalStorage<ord_t,val_t,Num,node_type> type;
    };

    //! Constructor
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    LocalStorage(const ordinal_type& sz,
		 const value_type& x = value_type(0.0)) { 
      c0 = x;
      c1 = x;
    }

    //! Default copy constructor

    //! Default destructor

    //! Default assignment operator

    //! Initialize values to a constant value
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_reference v) { 
      c0 = v;
      c1 = v;
    }

    //! Initialize values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      ordinal_type my_sz = sz;
      if (sz == 0) my_sz = Num;
      if (my_sz > 0) c0 = v[0];
      if (my_sz > 1) c1 = v[1];
    }

    //! Load values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void load(pointer v) { 
      v[0] = c0;
      v[1] = c1;
    }

    //! Resize to new size (values are preserved)
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void resize(const ordinal_type& sz) {}

    //! Reset storage to given array, size, and stride
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) {}

    //! Return size
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    static ordinal_type size() { return Num; }

    //! Coefficient access
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference operator[] (const ordinal_type& i) const { 
      if (i == 0) return c0;
      else if (i == 1) return c1;
      return c0;
    }

    //! Coefficient access
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference operator[] (const ordinal_type& i) { 
      if (i == 0) return c0;
      else if (i == 1) return c1;
      return c0;
    }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference getCoeff() {
      if (i == 0) return c0;
      else return c1;
    }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference getCoeff() const {
      if (i == 0) return c0;
      else return c1;
    }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_pointer coeff() const { return &c0; }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    pointer coeff() { return &c0; }

  private:

    //! Coefficient values
    value_type c0, c1;

  };

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t>
  class LocalStorage<ordinal_t, value_t, 4, KOKKOSARRAY_MACRO_DEVICE> {
  public:

    static const int Num = 4;
    
    static const bool is_static = true;
    static const int static_size = Num;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef KOKKOSARRAY_MACRO_DEVICE node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;

    //! Turn LocalStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef LocalStorage<ord_t,val_t,Num,node_type> type;
    };

    //! Constructor
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    LocalStorage(const ordinal_type& sz,
		 const value_type& x = value_type(0.0)) { 
      c0 = x;
      c1 = x;
      c2 = x;
      c3 = x;
    }

    //! Default copy constructor

    //! Default destructor

    //! Default assignment operator

    //! Initialize values to a constant value
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_reference v) { 
      c0 = v;
      c1 = v;
      c2 = v;
      c3 = v;
    }

    //! Initialize values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      ordinal_type my_sz = sz;
      if (sz == 0) my_sz = Num;
      if (my_sz > 0) c0 = v[0];
      if (my_sz > 1) c1 = v[1];
      if (my_sz > 2) c2 = v[2];
      if (my_sz > 3) c3 = v[3];
    }

    //! Load values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void load(pointer v) { 
      v[0] = c0;
      v[1] = c1;
      v[2] = c2;
      v[3] = c3;
    }

    //! Resize to new size (values are preserved)
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void resize(const ordinal_type& sz) {}

    //! Reset storage to given array, size, and stride
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) {}

    //! Return size
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    static ordinal_type size() { return Num; }

    //! Coefficient access
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference operator[] (const ordinal_type& i) const { 
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else if (i == 3) return c3;
      return c0;
    }

    //! Coefficient access
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference operator[] (const ordinal_type& i) { 
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else if (i == 3) return c3;
      return c0;
    }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference getCoeff() {
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else return c3;
    }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference getCoeff() const {
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else return c3;
    }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_pointer coeff() const { return &c0; }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    pointer coeff() { return &c0; }

  private:

    //! Coefficient values
    value_type c0, c1, c2, c3;

  };

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t>
  class LocalStorage<ordinal_t, value_t, 8, KOKKOSARRAY_MACRO_DEVICE> {
  public:

    static const int Num = 8;
    
    static const bool is_static = true;
    static const int static_size = Num;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef KOKKOSARRAY_MACRO_DEVICE node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;

    //! Turn LocalStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef LocalStorage<ord_t,val_t,Num,node_type> type;
    };

    //! Constructor
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    LocalStorage(const ordinal_type& sz,
		 const value_type& x = value_type(0.0)) { 
      c0 = x;
      c1 = x;
      c2 = x;
      c3 = x;
      c4 = x;
      c5 = x;
      c6 = x;
      c7 = x;
    }

    //! Default copy constructor

    //! Default destructor

    //! Default assignment operator

    //! Initialize values to a constant value
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_reference v) { 
      c0 = v;
      c1 = v;
      c2 = v;
      c3 = v;
      c4 = v;
      c5 = v;
      c6 = v;
      c7 = v;
    }

    //! Initialize values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      ordinal_type my_sz = sz;
      if (sz == 0) my_sz = Num;
      if (my_sz > 0) c0 = v[0];
      if (my_sz > 1) c1 = v[1];
      if (my_sz > 2) c2 = v[2];
      if (my_sz > 3) c3 = v[3];
      if (my_sz > 4) c4 = v[4];
      if (my_sz > 5) c5 = v[5];
      if (my_sz > 6) c6 = v[6];
      if (my_sz > 7) c7 = v[7];
    }

    //! Load values to an array of values
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void load(pointer v) { 
      v[0] = c0;
      v[1] = c1;
      v[2] = c2;
      v[3] = c3;
      v[4] = c4;
      v[5] = c5;
      v[6] = c6;
      v[7] = c7;
    }

    //! Resize to new size (values are preserved)
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void resize(const ordinal_type& sz) {}

    //! Reset storage to given array, size, and stride
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) {}

    //! Return size
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    static ordinal_type size() { return Num; }

    //! Coefficient access
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference operator[] (const ordinal_type& i) const { 
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else if (i == 3) return c3;
      else if (i == 4) return c4;
      else if (i == 5) return c5;
      else if (i == 6) return c6;
      else if (i == 7) return c7;
      return c0;
    }

    //! Coefficient access
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference operator[] (const ordinal_type& i) { 
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else if (i == 3) return c3;
      else if (i == 4) return c4;
      else if (i == 5) return c5;
      else if (i == 6) return c6;
      else if (i == 7) return c7;
      return c0;
    }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    reference getCoeff() {
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else if (i == 3) return c3;
      else if (i == 4) return c4;
      else if (i == 5) return c5;
      else if (i == 6) return c6;
      else return c7;
    }

    template <int i>
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference getCoeff() const {
      if (i == 0) return c0;
      else if (i == 1) return c1;
      else if (i == 2) return c2;
      else if (i == 3) return c3;
      else if (i == 4) return c4;
      else if (i == 5) return c5;
      else if (i == 6) return c6;
      else return c7;
    }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    const_pointer coeff() const { return &c0; }

    //! Get coefficients
    KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
    pointer coeff() { return &c0; }

  private:

    //! Coefficient values
    value_type c0, c1, c2, c3, c4, c5, c6, c7;

  };

}

#endif
