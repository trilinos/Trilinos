// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_CMATH_HPP
#define SACADO_CMATH_HPP

#include <cmath>        // for most math functions
#include "Sacado_ConfigDefs.h"

namespace Sacado {

  // Replacement for ternary operator, for scalar types that don't implement
  // logical operations that return bool, e.g., a simd scalar type that returns
  // a simd bool.  Sacado overloaded operators use this internally when ever
  // the ternary operator would be used.  It can also be used by client code.
  template <typename Cond, typename T>
  SACADO_INLINE_FUNCTION
  T if_then_else(const Cond cond, const T& a, const T& b) {
    return cond ? a : b;
  }

  // Special version of sqrt(x) that avoids the NaN if x==0 in the derivative.
  // The default implementation just calls the standard sqrt(x).
  template <typename T>
  SACADO_INLINE_FUNCTION
  T safe_sqrt(const T& x) {
    using std::sqrt;
    return sqrt(x);
  }

}

#endif // SACADO_CMATH_HPP
