// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_MP_EXPRESSIONTRAITS_HPP
#define SACADO_MP_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace MP {
    template <typename T> class Expr;
  }
}

namespace Sacado {

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< MP::Expr<T> > {
    typedef typename ScalarType< typename MP::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< MP::Expr<T> > {
    typedef typename MP::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< MP::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsScalarType< MP::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< MP::Expr<T> > {
    typedef typename ValueType< MP::Expr<T> >::type value_type;
    static const value_type& eval(const MP::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< MP::Expr<T> > {
    typedef typename ValueType< MP::Expr<T> >::type value_type;
    typedef typename ScalarType< MP::Expr<T> >::type scalar_type;
    static const scalar_type& eval(const MP::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_MP_EXPRESSIONTRAITS_HPP
