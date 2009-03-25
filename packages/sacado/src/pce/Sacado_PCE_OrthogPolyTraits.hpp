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

#ifndef SACADO_PCE_ORTHOGPOLYTRAITS_HPP
#define SACADO_PCE_ORTHOGPOLYTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace PCE {
    template <typename T> class OrthogPoly;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T>
  class Promote< PCE::OrthogPoly<T>, PCE::OrthogPoly<T> > {
  public:

    typedef PCE::OrthogPoly<T> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R>
  class Promote< PCE::OrthogPoly<L>, R > {
  public:

    typedef typename ValueType< PCE::OrthogPoly<L> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::OrthogPoly<value_type> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R>
  class Promote< L, PCE::OrthogPoly<R> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< PCE::OrthogPoly<R> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::OrthogPoly<value_type> type;
  };

  //! Specialization of %ScalarType to OrthogPoly types
  template <typename T>
  struct ScalarType< PCE::OrthogPoly<T> > {
    typedef typename ScalarType<typename PCE::OrthogPoly<T>::value_type>::type type;
  };

  //! Specialization of %ValueType to OrthogPoly types
  template <typename T>
  struct ValueType< PCE::OrthogPoly<T> > {
    typedef typename PCE::OrthogPoly<T>::value_type type;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T>
  struct IsADType< PCE::OrthogPoly<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T>
  struct IsScalarType< PCE::OrthogPoly<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to OrthogPoly types
  template <typename T>
  struct Value< PCE::OrthogPoly<T> > {
    typedef typename ValueType< PCE::OrthogPoly<T> >::type value_type;
    static const value_type& eval(const PCE::OrthogPoly<T>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to OrthogPoly types
  template <typename T>
  struct ScalarValue< PCE::OrthogPoly<T> > {
    typedef typename ValueType< PCE::OrthogPoly<T> >::type value_type;
    typedef typename ScalarType< PCE::OrthogPoly<T> >::type scalar_type;
    static const scalar_type& eval(const PCE::OrthogPoly<T>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to OrthogPoly types
  template <typename T>
  struct StringName< PCE::OrthogPoly<T> > {
    static std::string eval() { 
      return std::string("Sacado::PCE::OrthogPoly< ") + 
	StringName<T>::eval() + " >"; }
  };

} // namespace Sacado

#endif // SACADO_PCE_UNIVARIATEHERMITETRAITS_HPP
