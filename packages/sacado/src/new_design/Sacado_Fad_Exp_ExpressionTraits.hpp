// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_EXPRESSIONTRAITS_HPP
#define SACADO_FAD_EXP_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

#include "Sacado_Fad_Exp_Expression.hpp"

namespace Sacado {

  //! Specialization of %Promote to Expr types
  SACADO_EXPR_PROMOTE_SPEC( Fad::Exp )

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< Fad::Exp::Expr<T> > {
    typedef typename ScalarType< typename Fad::Exp::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< Fad::Exp::Expr<T> > {
    typedef typename Fad::Exp::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< Fad::Exp::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsScalarType to Expr types
  template <typename T>
  struct IsScalarType< Fad::Exp::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %IsSimdType to Expr types
  template <typename T>
  struct IsSimdType< Fad::Exp::Expr<T> > {
    static const bool value = IsSimdType< typename Fad::Exp::Expr<T>::value_type >::value;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< Fad::Exp::Expr<T> > {
    typedef typename ValueType< Fad::Exp::Expr<T> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const Fad::Exp::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< Fad::Exp::Expr<T> > {
    typedef typename ValueType< Fad::Exp::Expr<T> >::type value_type;
    typedef typename ScalarType< Fad::Exp::Expr<T> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const Fad::Exp::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_FAD_EXP_EXPRESSIONTRAITS_HPP
