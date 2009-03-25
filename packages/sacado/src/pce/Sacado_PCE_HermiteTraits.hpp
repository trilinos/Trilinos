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

#ifndef SACADO_PCE_HERMITETRAITS_HPP
#define SACADO_PCE_HERMITETRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace PCE {
    template <typename T> class Hermite;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T>
  class Promote< PCE::Hermite<T>, PCE::Hermite<T> > {
  public:

    typedef PCE::Hermite<T> type;
  };

  //! Specialization of %Promote to Hermite types
  template <typename L, typename R>
  class Promote< PCE::Hermite<L>, R > {
  public:

    typedef typename ValueType< PCE::Hermite<L> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::Hermite<value_type> type;
  };

  //! Specialization of %Promote to Hermite types
  template <typename L, typename R>
  class Promote< L, PCE::Hermite<R> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< PCE::Hermite<R> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::Hermite<value_type> type;
  };

  //! Specialization of %ScalarType to Hermite types
  template <typename T>
  struct ScalarType< PCE::Hermite<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ValueType to Hermite types
  template <typename T>
  struct ValueType< PCE::Hermite<T> > {
    typedef T type;
  };

  //! Specialization of %IsADType to Hermite types
  template <typename T>
  struct IsADType< PCE::Hermite<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Hermite types
  template <typename T>
  struct IsScalarType< PCE::Hermite<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Hermite types
  template <typename T>
  struct Value< PCE::Hermite<T> > {
    typedef typename ValueType< PCE::Hermite<T> >::type value_type;
    static const value_type& eval(const PCE::Hermite<T>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Hermite types
  template <typename T>
  struct ScalarValue< PCE::Hermite<T> > {
    typedef typename ValueType< PCE::Hermite<T> >::type value_type;
    typedef typename ScalarType< PCE::Hermite<T> >::type scalar_type;
    static const scalar_type& eval(const PCE::Hermite<T>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Hermite types
  template <typename T>
  struct StringName< PCE::Hermite<T> > {
    static std::string eval() { 
      return std::string("Sacado::PCE::Hermite< ") + 
	StringName<T>::eval() + " >"; }
  };

} // namespace Sacado

#endif // SACADO_PCE_HERMITETRAITS_HPP
