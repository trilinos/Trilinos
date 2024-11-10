// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_ETPCE_EXPRESSIONTRAITS_HPP
#define SACADO_ETPCE_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETPCE {
    template <typename T> class Expr;
  }
}

namespace Sacado {

  // We don't specialize Promote because otherwise we can get ambiguous
  // partial specializations with the OrthogPoly classes.

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< ETPCE::Expr<T> > {
    typedef typename ScalarType< typename ETPCE::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< ETPCE::Expr<T> > {
    typedef typename ETPCE::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< ETPCE::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsScalarType< ETPCE::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< ETPCE::Expr<T> > {
    typedef typename ValueType< ETPCE::Expr<T> >::type value_type;
    static const value_type& eval(const ETPCE::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< ETPCE::Expr<T> > {
    typedef typename ValueType< ETPCE::Expr<T> >::type value_type;
    typedef typename ScalarType< ETPCE::Expr<T> >::type scalar_type;
    static const scalar_type& eval(const ETPCE::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_ETPCE_EXPRESSIONTRAITS_HPP
