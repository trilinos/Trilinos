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

#ifndef SACADO_FAD_EXPRESSIONTRAITS_HPP
#define SACADO_FAD_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T> class Expr;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Expr types
  template <typename T>
  struct Promote< Fad::Expr<T>, Fad::Expr<T> > {
    typedef Fad::Expr<T> type;
  };

  //! Specialization of %Promote to Expr types
  template <typename T, typename R>
  struct Promote< Fad::Expr<T>, R > {
    typedef typename ValueType< Fad::Expr<T> >::type value_type_l;
    typedef typename Promote<R,R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::Expr<value_type> type;
  };

  //! Specialization of %Promote to Expr types
  template <typename L, typename T>
  struct Promote< L, Fad::Expr<T> > {
  public:

    typedef typename Promote<L,L>::type value_type_l;
    typedef typename ValueType< Fad::Expr<T> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::Expr<value_type> type;
  };

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< Fad::Expr<T> > {
    typedef T type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< Fad::Expr<T> > {
    typedef T type;
  };

   //! Specialization of %ScalarValueType to Expr types
  template <typename T>
  struct ScalarValueType< Fad::Expr<T> > {
    typedef typename ScalarValueType<T>::type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< Fad::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsScalarType< Fad::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< Fad::Expr<T> > {
    typedef typename ValueType< Fad::Expr<T> >::type value_type;
    static const value_type& eval(const Fad::Expr<T>& x) { 
      return x.val(); }
  };

} // namespace Sacado

#endif // SACADO_FAD_EXPRESSIONTRAITS_HPP
