// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PCE_UNIVARIATEHERMITETRAITS_HPP
#define SACADO_PCE_UNIVARIATEHERMITETRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace PCE {
    template <typename T> class UnivariateHermite;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T>
  class Promote< PCE::UnivariateHermite<T>, PCE::UnivariateHermite<T> > {
  public:

    typedef PCE::UnivariateHermite<T> type;
  };

  //! Specialization of %Promote to UnivariateHermite types
  template <typename L, typename R>
  class Promote< PCE::UnivariateHermite<L>, R > {
  public:

    typedef typename ValueType< PCE::UnivariateHermite<L> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::UnivariateHermite<value_type> type;
  };

  //! Specialization of %Promote to UnivariateHermite types
  template <typename L, typename R>
  class Promote< L, PCE::UnivariateHermite<R> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< PCE::UnivariateHermite<R> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::UnivariateHermite<value_type> type;
  };

  //! Specialization of %ScalarType to UnivariateHermite types
  template <typename T>
  struct ScalarType< PCE::UnivariateHermite<T> > {
    typedef T type;
  };

  //! Specialization of %ValueType to UnivariateHermite types
  template <typename T>
  struct ValueType< PCE::UnivariateHermite<T> > {
    typedef T type;
  };

   //! Specialization of %ScalarValueType to UnivariateHermite types
  template <typename T>
  struct ScalarValueType< PCE::UnivariateHermite<T> > {
    typedef typename ScalarValueType< T >::type type;
  };

  //! Specialization of %IsADType to UnivariateHermite types
  template <typename T>
  struct IsADType< PCE::UnivariateHermite<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to UnivariateHermite types
  template <typename T>
  struct IsScalarType< PCE::UnivariateHermite<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to UnivariateHermite types
  template <typename T>
  struct Value< PCE::UnivariateHermite<T> > {
    typedef typename ValueType< PCE::UnivariateHermite<T> >::type value_type;
    static const value_type& eval(const PCE::UnivariateHermite<T>& x) { 
      return x.val(); }
  };

} // namespace Sacado

#endif // SACADO_PCE_UNIVARIATEHERMITETRAITS_HPP
