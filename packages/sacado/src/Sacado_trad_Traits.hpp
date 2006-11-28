// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef SACADO_TRAD_TRAITS_HPP
#define SACADO_TRAD_TRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
#ifdef SACADO_NAMESPACE
#define SNS Sacado::Rad
namespace Sacado {
  namespace Rad {
    template <typename T> class ADvar;
    template <typename T> class ADvari;
  }
}
#else // SACADO_NAMESPACE  
#define SNS // nothing
template <typename T> class ADvar;
template <typename T> class ADvari;
#endif //SACADO_NAMESPACE  

namespace Sacado {

  //! Specialization of %Promote to ADvar types
  template <typename T>
  class Promote< SNS::ADvar<T>, SNS::ADvar<T> > {
  public:

    typedef SNS::ADvar<T> type;
  };

  //! Specialization of %Promote to ADvar types
  template <typename L, typename R>
  class Promote< SNS::ADvar<L>, R > {
  public:

    typedef typename ValueType< SNS::ADvar<L> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef SNS::ADvar<value_type> type;
  };

  //! Specialization of %Promote to ADvar types
  template <typename L, typename R>
  class Promote< L, SNS::ADvar<R> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< SNS::ADvar<R> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef SNS::ADvar<value_type> type;
  };

  //! Specialization of %ScalarType to ADvar types
  template <typename T>
  struct ScalarType< SNS::ADvar<T> > {
    typedef T type;
  };

  //! Specialization of %ScalarType to ADvari types
  template <typename T>
  struct ScalarType< SNS::ADvari<T> > {
    typedef T type;
  };

  //! Specialization of %ValueType to ADvar types
  template <typename T>
  struct ValueType< SNS::ADvar<T> > {
    typedef T type;
  };

  //! Specialization of %ValueType to ADvari types
  template <typename T>
  struct ValueType< SNS::ADvari<T> > {
    typedef T type;
  };

  //! Specialization of %ScalarValueType to ADvar types
  template <typename T>
  struct ScalarValueType< SNS::ADvar<T> > {
    typedef typename ScalarValueType< T >::type type;
  };

  //! Specialization of %ScalarValueType to ADvari types
  template <typename T>
  struct ScalarValueType< SNS::ADvari<T> > {
    typedef typename ScalarValueType< T >::type type;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsADType< SNS::ADvar<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsADType< SNS::ADvari<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsScalarType< SNS::ADvar<T> > {
    static const bool value = false;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsScalarType< SNS::ADvari<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to ADvar types
  template <typename T>
  struct Value< SNS::ADvar<T> > {
    typedef typename ValueType< SNS::ADvar<T> >::type value_type;
    static const value_type& eval(const SNS::ADvar<T>& x) { 
      return x.val(); }
  };

  //! Specialization of %Value to ADvari types
  template <typename T>
  struct Value< SNS::ADvari<T> > {
    typedef typename ValueType< SNS::ADvari<T> >::type value_type;
    static const value_type& eval(const SNS::ADvari<T>& x) { 
      return x.val(); }
  };

} // namespace Sacado

#undef SNS

#endif // SACADO_TRAD_TRAITS_HPP
