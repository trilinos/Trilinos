// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_LFAD_EXPRESSIONTRAITS_HPP
#define SACADO_LFAD_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace LFad {
    template <typename T> class Expr;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Expr types
  SACADO_EXPR_PROMOTE_SPEC( LFad )

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< LFad::Expr<T> > {
    typedef typename ScalarType< typename LFad::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< LFad::Expr<T> > {
    typedef typename LFad::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< LFad::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsScalarType< LFad::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< LFad::Expr<T> > {
    typedef typename ValueType< LFad::Expr<T> >::type value_type;
    static const value_type& eval(const LFad::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< LFad::Expr<T> > {
    typedef typename ValueType< LFad::Expr<T> >::type value_type;
    typedef typename ScalarType< LFad::Expr<T> >::type scalar_type;
    static const scalar_type& eval(const LFad::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_LFAD_EXPRESSIONTRAITS_HPP
