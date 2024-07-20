// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_ELRFAD_EXPRESSIONTRAITS_HPP
#define SACADO_ELRFAD_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ELRFad {
    template <typename T> class Expr;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Expr types
  SACADO_EXPR_PROMOTE_SPEC( ELRFad )

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< ELRFad::Expr<T> > {
    typedef typename ScalarType< typename ELRFad::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< ELRFad::Expr<T> > {
    typedef typename ELRFad::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< ELRFad::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsScalarType< ELRFad::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< ELRFad::Expr<T> > {
    typedef typename ValueType< ELRFad::Expr<T> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const ELRFad::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< ELRFad::Expr<T> > {
    typedef typename ValueType< ELRFad::Expr<T> >::type value_type;
    typedef typename ScalarType< ELRFad::Expr<T> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const ELRFad::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_ELRFAD_EXPRESSIONTRAITS_HPP
