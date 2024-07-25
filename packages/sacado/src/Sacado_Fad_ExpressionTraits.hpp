// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXPRESSIONTRAITS_HPP
#define SACADO_FAD_EXPRESSIONTRAITS_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_ExpressionTraits.hpp"

#else

#include "Sacado_Traits.hpp"

#include "Sacado_Fad_ExpressionFwd.hpp"

namespace Sacado {

  //! Specialization of %Promote to Expr types
  SACADO_EXPR_PROMOTE_SPEC( Fad )

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< Fad::Expr<T> > {
    typedef typename ScalarType< typename Fad::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< Fad::Expr<T> > {
    typedef typename Fad::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< Fad::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsSclarType to Expr types
  template <typename T>
  struct IsScalarType< Fad::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %IsSimdType to Expr types
  template <typename T>
  struct IsSimdType< Fad::Expr<T> > {
    static const bool value = IsSimdType< typename Fad::Expr<T>::value_type >::value;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< Fad::Expr<T> > {
    typedef typename ValueType< Fad::Expr<T> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const Fad::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< Fad::Expr<T> > {
    typedef typename ValueType< Fad::Expr<T> >::type value_type;
    typedef typename ScalarType< Fad::Expr<T> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const Fad::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_EXPRESSIONTRAITS_HPP
