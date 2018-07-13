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

#ifndef SACADO_CMATH_HPP
#define SACADO_CMATH_HPP

#include <cmath>        // for most math functions
#include "Sacado_ConfigDefs.h"

// Define some math functions that aren't usually in cmath
#if !( (defined(_GLIBCXX_USE_C99_MATH_TR1) && defined(__GXX_EXPERIMENTAL_CXX0X__)) || defined(HAVE_SACADO_CXX11) || defined(HAS_C99_TR1_CMATH) || defined(USER_DISABLE_SACADO_TR1_CMATH) )
namespace std {
  inline float acosh(float x) {
    return std::log(x + std::sqrt(x*x - float(1.0))); }
  inline float asinh(float x) {
    return std::log(x + std::sqrt(x*x + float(1.0))); }
  inline float atanh(float x) {
    return float(0.5)*std::log((float(1.0)+x)/(float(1.0)-x)); }

  inline double acosh(double x) {
    return std::log(x + std::sqrt(x*x - double(1.0))); }
  inline double asinh(double x) {
    return std::log(x + std::sqrt(x*x + double(1.0))); }
  inline double atanh(double x) {
    return double(0.5)*std::log((double(1.0)+x)/(double(1.0)-x)); }
}
#endif // HAS_C99_TR1_CMATH

namespace Sacado {

  // Replacement for ternary operator, for scalar types that don't implement
  // logical operations that return bool, e.g., a simd scalar type that returns
  // a simd bool.  Sacado overloaded operators use this internally when ever
  // the ternary operator would be used.  It can also be used by client code.
  template <typename Cond, typename T>
  KOKKOS_INLINE_FUNCTION
  T if_then_else(const Cond cond, const T& a, const T& b) {
    return cond ? a : b;
  }

}

#endif // SACADO_CMATH_HPP
