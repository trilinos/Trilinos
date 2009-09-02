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

#ifndef SACADO_TAY_TAYLORTRAITS_HPP
#define SACADO_TAY_TAYLORTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Tay {
    template <typename T> class Taylor;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SimpleTaylor types
  template <typename T>
  class Promote< Tay::Taylor<T>, Tay::Taylor<T> > {
  public:

    typedef Tay::Taylor<T> type;
  };

  //! Specialization of %Promote to SimpleTaylor types
  template <typename L, typename R>
  class Promote< Tay::Taylor<L>, R > {
  public:

    typedef typename ValueType< Tay::Taylor<L> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Tay::Taylor<value_type> type;
  };

  //! Specialization of %Promote to SimpleTaylor types
  template <typename L, typename R>
  class Promote< L, Tay::Taylor<R> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< Tay::Taylor<R> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Tay::Taylor<value_type> type;
  };

  //! Specialization of %ScalarType to Taylor types
  template <typename T>
  struct ScalarType< Tay::Taylor<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ValueType to Taylor types
  template <typename T>
  struct ValueType< Tay::Taylor<T> > {
    typedef T type;
  };

  //! Specialization of %IsADType to Taylor types
  template <typename T>
  struct IsADType< Tay::Taylor<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Taylor types
  template <typename T>
  struct IsScalarType< Tay::Taylor<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Taylor types
  template <typename T>
  struct Value< Tay::Taylor<T> > {
    typedef typename ValueType< Tay::Taylor<T> >::type value_type;
    static const value_type& eval(const Tay::Taylor<T>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Taylor types
  template <typename T>
  struct ScalarValue< Tay::Taylor<T> > {
    typedef typename ValueType< Tay::Taylor<T> >::type value_type;
    typedef typename ScalarType< Tay::Taylor<T> >::type scalar_type;
    static const scalar_type& eval(const Tay::Taylor<T>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Taylor types
  template <typename T>
  struct StringName< Tay::Taylor<T> > {
    static std::string eval() { 
      return std::string("Sacado::Tay::Taylor< ") + 
	StringName<T>::eval() + " >"; }
  };

} // namespace Sacado

#endif // SACADO_TAYLOR_SIMPLETAYLORTRAITS_HPP
