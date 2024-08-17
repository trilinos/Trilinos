//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_ARITHTRAITS_HPP
#define KOKKOS_ARITHTRAITS_HPP

/// \file Kokkos_ArithTraits.hpp
/// \brief Declaration and definition of Kokkos::ArithTraits

#include <KokkosKernels_config.h>
#include <Kokkos_NumericTraits.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_Macros.hpp>

#include <impl/Kokkos_QuadPrecisionMath.hpp>

#include <cfloat>
#include <climits>
#include <cmath>
#include <complex>  // std::complex
#include <limits>   // std::numeric_limits
#ifdef __CUDACC__
#include <math_constants.h>
#endif

namespace {  // anonymous

/// \fn intPowImpl
/// \tparam IntType A built-in integer type.
/// \brief Implementation of intPowSigned and intPowUnsigned.
///
/// \pre x != 0
/// \pre y > 0
///
/// Use intPowSigned or intPowUnsigned for general y.
template <class IntType>
KOKKOS_FORCEINLINE_FUNCTION IntType intPowImpl(const IntType x, const IntType y) {
  // Recursion (unrolled into while loop): pow(x, 2y) = (x^y)^2
  IntType prod  = x;
  IntType y_cur = 1;
  // If y == 1, then prod stays x.
  while (y_cur < y) {
    prod  = prod * prod;
    y_cur = y_cur << 1;
  }
  // abs(y - y_cur) < floor(log2(y)), so it won't hurt asymptotic run
  // time to finish the remainder in a linear iteration.
  if (y > y_cur) {
    const IntType left = y - y_cur;
    for (IntType k = 0; k < left; ++k) {
      prod = prod * x;
    }
  } else if (y < y_cur) {
    // There's probably a better way to do this in order to avoid the
    // (expensive) integer division, but I'm not motivated to think of
    // it at the moment.
    const IntType left = y_cur - y;
    for (IntType k = 0; k < left; ++k) {
      prod = prod / x;
    }
  }
  return prod;

  // y = 8:
  //
  // x,1   -> x^2,2
  // x^2,2 -> x^4,4
  // x^4,4 -> x^8,8
  //
  // y = 9:
  //
  // x,1   -> x^2,2
  // x^2,2 -> x^4,4
  // x^4,4 -> x^8,8
  //
  // y - y_cur is what's left over.  Just do it one at a time.
  //
  // y = 3:
  // x,1   -> x^2,2
  // x^2,2 -> x^4,4
}

// Warning free abs function for types where we don't know whether they are
// signed (like char)
template <class T, bool is_signed = std::numeric_limits<T>::is_signed>
struct integer_abs {
  static KOKKOS_INLINE_FUNCTION T abs(const T& val);
};

template <class T>
struct integer_abs<T, true> {
  static KOKKOS_INLINE_FUNCTION T abs(const T& x) { return x < 0 ? -x : x; }
};

template <class T>
struct integer_abs<T, false> {
  static KOKKOS_INLINE_FUNCTION T abs(const T& x) { return x; }
};

/// \fn intPowSigned
/// \tparam IntType A built-in signed integer type.
/// \brief Compute x raised to the power y.
///
/// If the arguments are invalid (e.g., if x and y are both zero), the
/// result of this function is undefined.  However, this function will
/// not throw an exception in that case.
template <class IntType>
KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<std::numeric_limits<IntType>::is_signed, IntType>::type
intPowSigned(const IntType x, const IntType y) {
  // It's not entirely clear what to return if x and y are both zero.
  // In the case of floating-point numbers, 0^0 is NaN.  Here, though,
  // I think it's safe to return 0.
  if (x == 0) {
    return 0;
  } else if (y == 0) {
    return 1;
  } else if (y < 0) {
    if (x == 1) {
      return 1;
    } else if (x == -1) {
      return (y % 2 == 0) ? 1 : -1;
    } else {
      return 0;  // round the fraction to zero
    }
  }
  return intPowImpl<IntType>(x, y);
}
template <class IntType>
KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<!std::numeric_limits<IntType>::is_signed, IntType>::type
intPowSigned(const IntType x, const IntType y) {
  // It's not entirely clear what to return if x and y are both zero.
  // In the case of floating-point numbers, 0^0 is NaN.  Here, though,
  // I think it's safe to return 0.
  if (x == 0) {
    return 0;
  } else if (y == 0) {
    return 1;
  }
  return intPowImpl<IntType>(x, y);
}

/// \fn intPowUnsigned
/// \tparam IntType A built-in unsigned integer type.
/// \brief Compute x raised to the power y.
///
/// If the arguments are invalid (e.g., if x and y are both zero), the
/// result of this function is undefined.  However, this function will
/// not throw an exception in that case.
template <class IntType>
KOKKOS_FORCEINLINE_FUNCTION IntType intPowUnsigned(const IntType x, const IntType y) {
  // It's not entirely clear what to return if x and y are both zero.
  // In the case of floating-point numbers, 0^0 is NaN.  Here, though,
  // I think it's safe to return 0.
  if (x == 0) {
    return 0;
  } else if (y == 0) {
    return 1;
  } else {
    return intPowImpl<IntType>(x, y);
  }
}

// It might make sense to use special sqrt() approximations for
// integer arguments, like those presented on the following web site:
//
// http://www.azillionmonkeys.com/qed/sqroot.html#implementations
//
// Note that some of the implementations on the above page break ANSI
// C(++) aliasing rules (by assigning to the results of
// reinterpret_cast-ing between int and float).  It's also just a
// performance optimization and not required for a reasonable
// implementation.

}  // namespace

namespace Kokkos {

// Macro to automate the wrapping of Kokkos Mathematical Functions
#define KOKKOSKERNELS_ARITHTRAITS_REAL_FP(FUNC_QUAL)                                               \
  static FUNC_QUAL val_type zero() { return static_cast<val_type>(0); }                            \
  static FUNC_QUAL val_type one() { return static_cast<val_type>(1); }                             \
  static FUNC_QUAL val_type min() { return Kokkos::Experimental::finite_min<val_type>::value; }    \
  static FUNC_QUAL val_type max() { return Kokkos::Experimental::finite_max<val_type>::value; }    \
  static FUNC_QUAL val_type infinity() { return Kokkos::Experimental::infinity<val_type>::value; } \
  static FUNC_QUAL val_type nan() { return Kokkos::Experimental::quiet_NaN<val_type>::value; }     \
  static FUNC_QUAL mag_type epsilon() { return Kokkos::Experimental::epsilon<val_type>::value; }   \
  static FUNC_QUAL mag_type sfmin() { return Kokkos::Experimental::norm_min<val_type>::value; }    \
  static FUNC_QUAL int base() { return Kokkos::Experimental::radix<val_type>::value; }             \
  static FUNC_QUAL mag_type prec() { return epsilon() * static_cast<mag_type>(base()); }           \
  static FUNC_QUAL int t() { return Kokkos::Experimental::digits<val_type>::value; }               \
  static FUNC_QUAL mag_type rnd() { return one(); }                                                \
  static FUNC_QUAL int emin() { return Kokkos::Experimental::min_exponent<val_type>::value; }      \
  static FUNC_QUAL mag_type rmin() { return Kokkos::Experimental::norm_min<val_type>::value; }     \
  static FUNC_QUAL int emax() { return Kokkos::Experimental::max_exponent<val_type>::value; }      \
  static FUNC_QUAL mag_type rmax() { return Kokkos::Experimental::finite_max<val_type>::value; }   \
                                                                                                   \
  static FUNC_QUAL bool isInf(const val_type x) { return Kokkos::isinf(x); }                       \
  static FUNC_QUAL bool isNan(const val_type x) { return Kokkos::isnan(x); }                       \
  static FUNC_QUAL mag_type abs(const val_type x) { return Kokkos::abs(x); }                       \
  static FUNC_QUAL mag_type real(const val_type x) { return Kokkos::real(x); }                     \
  static FUNC_QUAL mag_type imag(const val_type x) { return Kokkos::imag(x); }                     \
  static FUNC_QUAL val_type conj(const val_type x) { return x; }                                   \
  static FUNC_QUAL val_type pow(const val_type x, const val_type y) { return Kokkos::pow(x, y); }  \
  static FUNC_QUAL val_type sqrt(const val_type x) { return Kokkos::sqrt(x); }                     \
  static FUNC_QUAL val_type cbrt(const val_type x) { return Kokkos::cbrt(x); }                     \
  static FUNC_QUAL val_type exp(const val_type x) { return Kokkos::exp(x); }                       \
  static FUNC_QUAL val_type log(const val_type x) { return Kokkos::log(x); }                       \
  static FUNC_QUAL val_type log10(const val_type x) { return Kokkos::log10(x); }                   \
  static FUNC_QUAL val_type sin(const val_type x) { return Kokkos::sin(x); }                       \
  static FUNC_QUAL val_type cos(const val_type x) { return Kokkos::cos(x); }                       \
  static FUNC_QUAL val_type tan(const val_type x) { return Kokkos::tan(x); }                       \
  static FUNC_QUAL val_type sinh(const val_type x) { return Kokkos::sinh(x); }                     \
  static FUNC_QUAL val_type cosh(const val_type x) { return Kokkos::cosh(x); }                     \
  static FUNC_QUAL val_type tanh(const val_type x) { return Kokkos::tanh(x); }                     \
  static FUNC_QUAL val_type asin(const val_type x) { return Kokkos::asin(x); }                     \
  static FUNC_QUAL val_type acos(const val_type x) { return Kokkos::acos(x); }                     \
  static FUNC_QUAL val_type atan(const val_type x) { return Kokkos::atan(x); }                     \
                                                                                                   \
  static FUNC_QUAL bool isnaninf(const val_type x) { return isNan(x) || isInf(x); }                \
  static FUNC_QUAL magnitudeType magnitude(const val_type x) { return abs(x); }                    \
  static FUNC_QUAL val_type conjugate(const val_type x) { return conj(x); }                        \
  static FUNC_QUAL val_type squareroot(const val_type x) { return sqrt(x); }                       \
  static FUNC_QUAL mag_type eps() { return epsilon(); }

// Macro to automate the wrapping of Kokkos Mathematical Functions
#define KOKKOSKERNELS_ARITHTRAITS_HALF_FP(FUNC_QUAL)                                               \
  static FUNC_QUAL val_type zero() { return static_cast<val_type>(0); }                            \
  static FUNC_QUAL val_type one() { return static_cast<val_type>(1); }                             \
  static FUNC_QUAL val_type min() { return Kokkos::Experimental::finite_min<val_type>::value; }    \
  static FUNC_QUAL val_type max() { return Kokkos::Experimental::finite_max<val_type>::value; }    \
  static FUNC_QUAL val_type infinity() { return Kokkos::Experimental::infinity<val_type>::value; } \
  static FUNC_QUAL val_type nan() { return Kokkos::Experimental::quiet_NaN<val_type>::value; }     \
  static FUNC_QUAL mag_type epsilon() { return Kokkos::Experimental::epsilon<val_type>::value; }   \
  static FUNC_QUAL mag_type sfmin() { return Kokkos::Experimental::norm_min<val_type>::value; }    \
  static FUNC_QUAL int base() { return Kokkos::Experimental::radix<val_type>::value; }             \
  static FUNC_QUAL mag_type prec() { return epsilon() * static_cast<mag_type>(base()); }           \
  static FUNC_QUAL int t() { return Kokkos::Experimental::digits<val_type>::value; }               \
  static FUNC_QUAL mag_type rnd() { return one(); }                                                \
  static FUNC_QUAL int emin() { return Kokkos::Experimental::min_exponent<val_type>::value; }      \
  static FUNC_QUAL mag_type rmin() { return Kokkos::Experimental::norm_min<val_type>::value; }     \
  static FUNC_QUAL int emax() { return Kokkos::Experimental::max_exponent<val_type>::value; }      \
  static FUNC_QUAL mag_type rmax() { return Kokkos::Experimental::finite_max<val_type>::value; }   \
                                                                                                   \
  static FUNC_QUAL bool isInf(const val_type x) { return Kokkos::isinf(x); }                       \
  static FUNC_QUAL mag_type abs(const val_type x) { return Kokkos::abs(x); }                       \
  static FUNC_QUAL mag_type real(const val_type x) { return Kokkos::real(x); }                     \
  static FUNC_QUAL mag_type imag(const val_type x) { return Kokkos::imag(x); }                     \
  static FUNC_QUAL val_type conj(const val_type x) { return x; }                                   \
  static FUNC_QUAL val_type pow(const val_type x, const val_type y) { return Kokkos::pow(x, y); }  \
  static FUNC_QUAL val_type sqrt(const val_type x) { return Kokkos::sqrt(x); }                     \
  static FUNC_QUAL val_type cbrt(const val_type x) { return Kokkos::cbrt(x); }                     \
  static FUNC_QUAL val_type exp(const val_type x) { return Kokkos::exp(x); }                       \
  static FUNC_QUAL val_type log(const val_type x) { return Kokkos::log(x); }                       \
  static FUNC_QUAL val_type log10(const val_type x) { return Kokkos::log10(x); }                   \
  static FUNC_QUAL val_type sin(const val_type x) { return Kokkos::sin(x); }                       \
  static FUNC_QUAL val_type cos(const val_type x) { return Kokkos::cos(x); }                       \
  static FUNC_QUAL val_type tan(const val_type x) { return Kokkos::tan(x); }                       \
  static FUNC_QUAL val_type sinh(const val_type x) { return Kokkos::sinh(x); }                     \
  static FUNC_QUAL val_type cosh(const val_type x) { return Kokkos::cosh(x); }                     \
  static FUNC_QUAL val_type tanh(const val_type x) { return Kokkos::tanh(x); }                     \
  static FUNC_QUAL val_type asin(const val_type x) { return Kokkos::asin(x); }                     \
  static FUNC_QUAL val_type acos(const val_type x) { return Kokkos::acos(x); }                     \
  static FUNC_QUAL val_type atan(const val_type x) { return Kokkos::atan(x); }                     \
                                                                                                   \
  static FUNC_QUAL magnitudeType magnitude(const val_type x) { return abs(x); }                    \
  static FUNC_QUAL val_type conjugate(const val_type x) { return conj(x); }                        \
  static FUNC_QUAL val_type squareroot(const val_type x) { return sqrt(x); }                       \
  static FUNC_QUAL mag_type eps() { return epsilon(); }

#define KOKKOSKERNELS_ARITHTRAITS_CMPLX_FP(FUNC_QUAL)                                                                 \
                                                                                                                      \
  static constexpr bool is_specialized = true;                                                                        \
  static constexpr bool is_signed      = true;                                                                        \
  static constexpr bool is_integer     = false;                                                                       \
  static constexpr bool is_exact       = false;                                                                       \
  static constexpr bool is_complex     = true;                                                                        \
  static constexpr bool has_infinity   = true;                                                                        \
                                                                                                                      \
  using magnitudeType   = mag_type;                                                                                   \
  using halfPrecision   = ::Kokkos::complex<ArithTraits<mag_type>::halfPrecision>;                                    \
  using doublePrecision = ::Kokkos::complex<ArithTraits<mag_type>::doublePrecision>;                                  \
                                                                                                                      \
  static constexpr bool isComplex            = true;                                                                  \
  static constexpr bool isOrdinal            = false;                                                                 \
  static constexpr bool isComparable         = false;                                                                 \
  static constexpr bool hasMachineParameters = ArithTraits<mag_type>::hasMachineParameters;                           \
                                                                                                                      \
  static FUNC_QUAL val_type zero() { return val_type(ArithTraits<mag_type>::zero(), ArithTraits<mag_type>::zero()); } \
  static FUNC_QUAL val_type one() { return val_type(ArithTraits<mag_type>::one(), ArithTraits<mag_type>::zero()); }   \
  static FUNC_QUAL val_type min() { return val_type(ArithTraits<mag_type>::min(), ArithTraits<mag_type>::min()); }    \
  static FUNC_QUAL val_type max() { return val_type(ArithTraits<mag_type>::max(), ArithTraits<mag_type>::max()); }    \
  static FUNC_QUAL val_type infinity() {                                                                              \
    return val_type(ArithTraits<mag_type>::infinity(), ArithTraits<mag_type>::infinity());                            \
  }                                                                                                                   \
  static FUNC_QUAL val_type nan() { return val_type(ArithTraits<mag_type>::nan(), ArithTraits<mag_type>::nan()); }    \
  static FUNC_QUAL mag_type epsilon() { return ArithTraits<mag_type>::epsilon(); }                                    \
  static FUNC_QUAL mag_type sfmin() { return ArithTraits<mag_type>::sfmin(); }                                        \
  static FUNC_QUAL int base() { return ArithTraits<mag_type>::base(); }                                               \
  static FUNC_QUAL mag_type prec() { return ArithTraits<mag_type>::prec(); }                                          \
  static FUNC_QUAL int t() { return ArithTraits<mag_type>::t(); }                                                     \
  static FUNC_QUAL mag_type rnd() { return ArithTraits<mag_type>::rnd(); }                                            \
  static FUNC_QUAL int emin() { return ArithTraits<mag_type>::emin(); }                                               \
  static FUNC_QUAL mag_type rmin() { return ArithTraits<mag_type>::rmin(); }                                          \
  static FUNC_QUAL int emax() { return ArithTraits<mag_type>::emax(); }                                               \
  static FUNC_QUAL mag_type rmax() { return ArithTraits<mag_type>::rmax(); }                                          \
  static FUNC_QUAL bool isInf(const val_type x) {                                                                     \
    return ArithTraits<mag_type>::isInf(x.real()) || ArithTraits<mag_type>::isInf(x.imag());                          \
  }                                                                                                                   \
  static FUNC_QUAL bool isNan(const val_type x) {                                                                     \
    return ArithTraits<mag_type>::isNan(x.real()) || ArithTraits<mag_type>::isNan(x.imag());                          \
  }                                                                                                                   \
  static FUNC_QUAL mag_type abs(const val_type x) { return ::Kokkos::abs(x); }                                        \
  static FUNC_QUAL mag_type real(const val_type x) { return x.real(); }                                               \
  static FUNC_QUAL mag_type imag(const val_type x) { return x.imag(); }                                               \
  static FUNC_QUAL val_type conj(const val_type x) { return ::Kokkos::conj(x); }                                      \
  static FUNC_QUAL val_type pow(const val_type x, const val_type y) { return Kokkos::pow(x, y); }                     \
  static FUNC_QUAL val_type pow(const val_type x, const mag_type y) { return Kokkos::pow(x, y); }                     \
  static FUNC_QUAL val_type pow(const mag_type x, const val_type y) { return Kokkos::pow(x, y); }                     \
  static FUNC_QUAL val_type sqrt(const val_type x) { return ::Kokkos::sqrt(x); }                                      \
  static FUNC_QUAL val_type exp(const val_type x) { return Kokkos::exp(x); }                                          \
  static FUNC_QUAL val_type log(const val_type x) { return Kokkos::log(x); }                                          \
  static FUNC_QUAL val_type log10(const val_type x) { return Kokkos::log10(x); }                                      \
  static FUNC_QUAL val_type sin(const val_type x) { return Kokkos::sin(x); }                                          \
  static FUNC_QUAL val_type cos(const val_type x) { return Kokkos::cos(x); }                                          \
  static FUNC_QUAL val_type tan(const val_type x) { return Kokkos::tan(x); }                                          \
  static FUNC_QUAL val_type sinh(const val_type x) { return Kokkos::sinh(x); }                                        \
  static FUNC_QUAL val_type cosh(const val_type x) { return Kokkos::cosh(x); }                                        \
  static FUNC_QUAL val_type tanh(const val_type x) { return Kokkos::tanh(x); }                                        \
  static FUNC_QUAL val_type asin(const val_type x) { return Kokkos::asin(x); }                                        \
  static FUNC_QUAL val_type acos(const val_type x) { return Kokkos::acos(x); }                                        \
  static FUNC_QUAL val_type atan(const val_type x) { return Kokkos::atan(x); }                                        \
  static FUNC_QUAL bool isnaninf(const val_type& x) { return isNan(x) || isInf(x); }                                  \
  static FUNC_QUAL mag_type magnitude(const val_type x) { return abs(x); }                                            \
  static FUNC_QUAL val_type conjugate(const val_type x) { return conj(x); }                                           \
  static FUNC_QUAL val_type squareroot(const val_type x) { return sqrt(x); }                                          \
  static FUNC_QUAL mag_type eps() { return epsilon(); }

template <typename val_type>
static KOKKOS_FUNCTION typename std::enable_if<std::numeric_limits<val_type>::is_signed, val_type>::type
KokkosKernelsAbs(const val_type x) {
  return Kokkos::abs(x);
}

template <typename val_type>
static KOKKOS_FUNCTION typename std::enable_if<!std::numeric_limits<val_type>::is_signed, val_type>::type
KokkosKernelsAbs(const val_type x) {
  return x;
}

template <typename val_type>
static KOKKOS_FUNCTION typename std::enable_if<std::numeric_limits<val_type>::is_signed, val_type>::type
KokkosKernelsNan() {
  return -1;
}

template <typename val_type>
static KOKKOS_FUNCTION typename std::enable_if<!std::numeric_limits<val_type>::is_signed, val_type>::type
KokkosKernelsNan() {
  return Kokkos::Experimental::finite_max<val_type>::value;
}

#define KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()                                                                       \
                                                                                                                   \
  static constexpr bool is_specialized = true;                                                                     \
  static constexpr bool is_integer     = true;                                                                     \
  static constexpr bool is_exact       = true;                                                                     \
  static constexpr bool is_complex     = false;                                                                    \
  static constexpr bool has_infinity   = false;                                                                    \
                                                                                                                   \
  using magnitudeType   = mag_type;                                                                                \
  using halfPrecision   = val_type;                                                                                \
  using doublePrecision = val_type;                                                                                \
                                                                                                                   \
  static constexpr bool isComplex            = false;                                                              \
  static constexpr bool isOrdinal            = true;                                                               \
  static constexpr bool isComparable         = true;                                                               \
  static constexpr bool hasMachineParameters = false;                                                              \
                                                                                                                   \
  static KOKKOS_FUNCTION val_type zero() { return static_cast<val_type>(0); }                                      \
  static KOKKOS_FUNCTION val_type one() { return static_cast<val_type>(1); }                                       \
  static KOKKOS_FUNCTION val_type min() { return Kokkos::Experimental::finite_min<val_type>::value; }              \
  static KOKKOS_FUNCTION val_type max() { return Kokkos::Experimental::finite_max<val_type>::value; }              \
  static KOKKOS_FUNCTION val_type infinity() { return static_cast<val_type>(0); }                                  \
  static KOKKOS_FUNCTION val_type nan() { return KokkosKernelsNan<val_type>(); }                                   \
  static KOKKOS_FUNCTION bool isInf(const val_type) { return false; }                                              \
  static KOKKOS_FUNCTION bool isNan(const val_type) { return false; }                                              \
  static KOKKOS_FUNCTION mag_type abs(const val_type x) { return KokkosKernelsAbs(x); }                            \
  static KOKKOS_FUNCTION mag_type real(const val_type x) { return Kokkos::real(x); }                               \
  static KOKKOS_FUNCTION mag_type imag(const val_type) { return zero(); }                                          \
  static KOKKOS_FUNCTION val_type conj(const val_type x) { return x; }                                             \
  static KOKKOS_FUNCTION val_type pow(const val_type x, const val_type y) { return Kokkos::pow(x, y); }            \
  static KOKKOS_FUNCTION val_type sqrt(const val_type x) { return static_cast<val_type>(Kokkos::sqrt(abs(x))); }   \
  static KOKKOS_FUNCTION val_type cbrt(const val_type x) { return static_cast<val_type>(Kokkos::cbrt(abs(x))); }   \
  static KOKKOS_FUNCTION val_type exp(const val_type x) { return static_cast<val_type>(Kokkos::exp(abs(x))); }     \
  static KOKKOS_FUNCTION val_type log(const val_type x) { return static_cast<val_type>(Kokkos::log(abs(x))); }     \
  static KOKKOS_FUNCTION val_type log10(const val_type x) { return static_cast<val_type>(Kokkos::log10(abs(x))); } \
  static KOKKOS_FUNCTION mag_type epsilon() { return zero(); }                                                     \
  static KOKKOS_FUNCTION magnitudeType magnitude(const val_type x) { return abs(x); }                              \
  static KOKKOS_FUNCTION val_type conjugate(const val_type x) { return conj(x); }                                  \
  static KOKKOS_FUNCTION bool isnaninf(const val_type) { return false; }                                           \
  static KOKKOS_FUNCTION val_type squareroot(const val_type x) { return sqrt(x); }

/// \class ArithTraits
/// \brief Traits class for arithmetic on type T.
/// \tparam T "Scalar" type of interest
///
/// This is a traits class for the "arithmetic" type T.  "Arithmetic
/// types" include built-in signed and unsigned integer types,
/// floating-point types, complex-valued types, and anything else that
/// looks like these.  This class is useful for implementing numerical
/// algorithms that are generic on the data type.  You may also use
/// this class to query attributes of T, like whether it is signed or
/// complex, or its precision.
///
/// We really did not want to implement this class or expose it to
/// users.  It would be much better to use existing traits classes
/// like std::numeric_limits.  We decided to implement and expose this
/// class for the following reasons:
/// <ol>
/// <li> std::numeric_limits class methods cannot be used in CUDA
///      device functions, since they themselves are not device
///      functions </li>
/// <li> Existing traits classes like std::numeric_limits do not
///      provide enough information to implement algorithms that are
///      agnostic of whether T is real-valued or complex-valued. </li>
/// </ol>
///
/// All class methods must be suitable for parallel kernels, if the
/// type T itself is suitable for parallel kernels.  In particular,
/// specializations for types T that make sense to use on a CUDA
/// device must mark all class methods as device (and host) functions,
/// using the KOKKOS_FORCEINLINE_FUNCTION macro.  All class methods must be
/// callable both inside and outside a parallel kernel (for CUDA, this
/// means they must be marked as both device and host functions).
///
/// \section Kokkos_ArithTraits_compat Compatibility
///
/// Whenever possible, class methods in ArithTraits use the same names
/// as their equivalents in the C++ Standard Library.  If this was not
/// possible, for example with isInf and isNan, we explain why in
/// their documentation.
///
/// This class has redundant typedefs and methods in order to maintain
/// backwards compatibility with Teuchos::ScalarTraits, while
/// preferring forwards (partial) compatibility with
/// std::numeric_limits.  Users should prefer typedefs, \c bool
/// constants, and class methods compatible with std::numeric_limits,
/// to those from Teuchos::ScalarTraits.  The latter may go away at
/// any time.  Furthermore, Teuchos::ScalarTraits contains methods
/// that do not make sense for use as parallel device functions, in
/// particular those relating to pseudorandom number generation that
/// refer to hidden state, so we will never include all class methods
/// from Teuchos::ScalarTraits in ArithTraits.
///
/// \section Kokkos_ArithTraits_unsupp Unsupported types on CUDA devices
///
/// CUDA does not support long double or std::complex<T> in device
/// functions.  ArithTraits does have specializations for these types,
/// but the class methods therein are not marked as device functions.
///
/// \section Kokkos_ArithTraits_whyNotC99 What about C99 integer types?
///
/// C99 and C++11 include typedefs int${N}_t and uint${N}_t, where N
/// is the number of bits in the integer.  These typedefs are useful
/// because they make the length of the type explicit.  Users are
/// welcome to use these types as the template parameter of
/// ArithTraits.
///
/// We chose not to use these types when <i>defining</i> full
/// specializations of ArithTraits.  This is because the C99 integer
/// types are typedefs, not types in themselves.  This makes it
/// impossible to avoid duplicate or missing full specializations of
/// ArithTraits.  For example, on my Mac, for CUDA 5.5, gcc 4.2.1, and
/// Clang 3.2, <tt>int64_t</tt> is a typedef of <tt>long long</tt>,
/// but <tt>long long</tt> and <tt>long</tt> are separate types, even
/// though they have the same length (64 bits).  In contrast, on
/// Windows (even Win64), <tt>long</tt> is a 32-bit type (but a
/// distinct type from <tt>int</tt>), and <tt>long long</tt> is a
/// 64-bit type.  Thus, if we define full specializations of
/// ArithTraits using <i>only</i> the C99 integer types, we will be
/// missing a specialization for <tt>long</tt> on at least one
/// platform.
///
/// Rather than trouble ourselves with trying to figure this out for
/// each platform, we decided to provide specializations only for the
/// integer types in the C89 and C++03 language standards.  This
/// includes signed and unsigned versions of <tt>char</tt>,
/// <tt>short</tt>, <tt>int</tt>, and <tt>long</tt>.  We also include
/// <tt>long long</tt> if your platform supports it.  We may thus have
/// left out some C99 integer type, but this is only possible if the
/// C89 / C++03 integer types do not have complete coverage of all
/// powers of two bits from 8 up to the longest provided length (e.g.,
/// 64 on a 64-bit system).  On all platforms I have encountered,
/// <tt>char</tt> has 8 bits and <tt>short</tt> has 16 bits, so I am
/// not worried about missing specializations for <tt>int16_t</tt> or
/// <tt>uint16_t</tt>.  If you should find that either of these
/// specializations are missing, though, please let us know.
///
/// Note that <tt>char</tt>, <tt>signed char</tt>, and <tt>unsigned
/// char</tt> are distinct types, whether <tt>char</tt> is signed or
/// unsigned.  (The language standards do not specify whether
/// <tt>char</tt> is signed or unsigned.)  That is, <tt>char</tt> is
/// <i>not</i> a typedef of <tt>signed char</tt> or <tt>unsigned
/// char</tt>.  This is why we provide full specializations of
/// ArithTraits for each of these types.  Interestingly enough, on my
/// system, <tt>char</tt> and <tt>int8_t</tt> are different types, but
/// <tt>signed char</tt> and <tt>int8_t</tt> are the same.
///
/// \section Kokkos_ArithTraits_impl Implementation notes
///
/// This section contains notes to developers who which to add a
/// partial specialization of this class for a new type T.  If you
/// decide to write a default templated implementation, it must not
/// declare any methods as device functions.  This ensures correct
/// behavior for arbitrary T, but does require specializations for
/// common types like T = float and double, as well as for other types
/// T that make sense to use on a CUDA device.
template <class T>
class ArithTraits {
 public:
  /// \brief A type that acts like T and works with Kokkos.
  ///
  /// This is usually just an alias for T.  However, some types T do
  /// not work well with Kokkos.  In that case, we use a mostly
  /// equivalent type here.  For example, ArithTraits<std::complex<R>
  /// >::val_type is Kokkos::complex<R>.
  using val_type = T;
  /// \brief The type of the magnitude (absolute value) of T.
  ///
  /// We define this as the type returned by abs() in this class.  If
  /// T is real (not complex), then \c val_type and \c mag_type are
  /// usually the same.  If T is <tt>std::complex<R></tt> for some R,
  /// then R and \c mag_type are usually the same.
  using mag_type = T;

  //! Whether ArithTraits has a specialization for T.
  static constexpr bool is_specialized = false;
  //! Whether T is a signed type (has negative values).
  static constexpr bool is_signed = false;
  //! Whether T is an integer type.
  static constexpr bool is_integer = false;
  /// \brief Whether T "uses exact representations."
  ///
  /// The opposite of is_exact is "is approximate," that is, "may
  /// commit rounding error."
  static constexpr bool is_exact = false;
  //! Whether T is a complex-valued type.
  static constexpr bool is_complex = false;

  /// \brief Whether x is Inf.
  ///
  /// This can only be true for floating-point types T that support
  /// Inf.  If T is a complex type, we say that a T instance x is Inf
  /// if and only if <tt>isinf(real(x)) || isinf(imag(x))</tt>.
  ///
  /// Unfortunately we can't call this "isinf" (the equivalent C99
  /// function), because CUDA appears to implement that function using
  /// a macro, rather than using a function (as C++11 requires).
  static KOKKOS_FUNCTION bool isInf(const T& x);

  /// \brief Whether x is NaN (not a number).
  ///
  /// This can only be true for floating-point types T that support
  /// NaN.  If T is a complex type, we say that a T instance x is NaN
  /// if and only if <tt>isNan(real(x)) || isNan(imag(x))</tt>.
  ///
  /// Unfortunately we can't call this "isnan" (the equivalent C99
  /// function), because CUDA appears to implement that function using
  /// a macro, rather than using a function (as C++11 requires).
  static KOKKOS_FUNCTION bool isNan(const T& x);

  //! The absolute value (magnitude) of x.
  static KOKKOS_FUNCTION mag_type abs(const T& x);

  //! The zero value of T; the arithmetic identity.
  static KOKKOS_FUNCTION T zero();

  //! The one value of T; the multiplicative identity.
  static KOKKOS_FUNCTION T one();

  /// \brief True if this type T is capable of representing the
  /// positive infinity as a distinct special value, as with
  /// std::numeric_limits<T>::has_infinity.
  static constexpr bool has_infinity = false;

  /// \brief Returns the special value "positive infinity", as
  /// represented by the floating-point type T. Only meaningful if
  /// KokkosArithTraits<T>::has_infinity == true. Provides same
  /// functionality as std::numeric_limits<T>::infinity().
  ///
  /// \note Would have liked to mark it as constexpr but then would
  /// not be able to provide the specialization for std::complex<T>
  /// since its constructor only becomes constexpr with C++14.
  static KOKKOS_FUNCTION T infinity();

  /// \brief The minimum possible value of T.
  ///
  /// If T is a real floating-point type, then this is the minimum
  /// <i>positive</i> value, as with std::numeric_limits<T>::min().
  static KOKKOS_FUNCTION T min();

  //! The maximum possible value of T.
  static KOKKOS_FUNCTION T max();

  /// \brief The real part of x.
  ///
  /// If \c is_complex is false, then this just returns x.
  static KOKKOS_FUNCTION mag_type real(const T& x);

  /// \brief The imaginary part of x.
  ///
  /// If \c is_complex is false, then this just returns zero().
  static KOKKOS_FUNCTION mag_type imag(const T&);

  /// \brief The complex conjugate of x.
  ///
  /// If \c is_complex is false, then this just returns x.
  static KOKKOS_FUNCTION T conj(const T&);

  //! x raised to the power y.
  static KOKKOS_FUNCTION T pow(const T& x, const T& y);

  /// \brief The square root of x.
  ///
  /// If T is an integer type, this is the floor of the square root.
  /// If T is a complex-valued type, then this method returns the
  /// principal branch of the square root.
  ///
  /// If T is real-valued and x is negative, the result of the square
  /// root is undefined in general.  (CUDA does not allow throwing
  /// exceptions in device functions.)  Implementations should return
  /// NaN if the type T supports this.  Of course, in that case, the
  /// square of the result will not equal x.
  static KOKKOS_FUNCTION T sqrt(const T& x);

  /// \brief The cubic root of x.
  ///
  /// If T is an integer type, this is the floor of the cubic root.
  /// If T is a complex-valued type, then this method returns the
  /// principal branch of the cubic root.
  ///
  /// If T is real-valued and x is negative, the result of the cubic
  /// root is undefined in general.  (CUDA does not allow throwing
  /// exceptions in device functions.)  Implementations should return
  /// NaN if the type T supports this.  Of course, in that case, the
  /// cubic of the result will not equal x.
  static KOKKOS_FUNCTION T cbrt(const T& x);

  /// \brief The natural (base e) exponential function of x.
  ///
  /// If T is an integer type, this is the floor of the exponential
  /// function.  If T is a complex-valued type, then this method
  /// returns \f$e^{x+iy} = e^x ( cos(y) + i sin(y) )\f$.
  ///
  static KOKKOS_FUNCTION T exp(const T& x);

  /// \brief The natural (base e) logarithm of x.
  ///
  /// If T is an integer type, this is the floor of the logarithm.  If
  /// T is a complex-valued type, then this method returns the
  /// principal branch of the logarithm.
  ///
  /// If T is real-valued and x is negative, the result of the
  /// logarithm is undefined in general.  (CUDA does not allow
  /// throwing exceptions in device functions.)  Implementations
  /// should return NaN if the type T supports this.  Of course, in
  /// that case, if y is the result, \f$e^y\f$ will not equal x.
  static KOKKOS_FUNCTION T log(const T& x);

  /// \brief The base ten logarithm of the input.
  ///
  /// If T is an integer type, this is the floor of the logarithm.  If
  /// T is a complex-valued type, then this method returns the
  /// principal branch of the logarithm.
  ///
  /// If T is real-valued and x is negative, the result of the
  /// logarithm is undefined in general.  (CUDA does not allow
  /// throwing exceptions in device functions.)  Implementations
  /// should return NaN if the type T supports this.  Of course, in
  /// that case, if y is the result, \f$10^y\f$ will not equal x.
  static KOKKOS_FUNCTION T log10(const T& x);

  /// Trigonometric and hyperbolic functions are not available
  /// for integer types. This is because asin(sin(x)) is not x
  /// when x is integer with a rounding error.
  ///
  ///  KJ: log, exp also has this problem. We probably need to
  ///      disable them for integer types instead of providing
  ///      functionality with floor.

  /// \brief The sin function of x
  ///
  static KOKKOS_FUNCTION T sin(const T& x);

  /// \brief The cos function of x
  ///
  static KOKKOS_FUNCTION T cos(const T& x);

  /// \brief The tan function of x
  ///
  static KOKKOS_FUNCTION T tan(const T& x);

  /// \brief The sin hyperbolic function of x
  ///
  static KOKKOS_FUNCTION T sinh(const T& x);

  /// \brief The cos hyperbolic function of x
  ///
  static KOKKOS_FUNCTION T cosh(const T& x);

  /// \brief The tan hyperbolic function of x
  ///
  static KOKKOS_FUNCTION T tanh(const T& x);

  /// \brief The asin function of x
  ///
  static KOKKOS_FUNCTION T asin(const T& x);

  /// \brief The acos function of x
  ///
  static KOKKOS_FUNCTION T acos(const T& x);

  /// \brief The atan function of x
  ///
  static KOKKOS_FUNCTION T atan(const T& x);

  /// \brief Return a silent NaN, if appropriate for T.
  ///
  /// If T does <i>not</i> implement a silent NaN, the return value is
  /// undefined, but calling this method is still allowed.
  static KOKKOS_FUNCTION T nan();

  /// \brief Machine epsilon.
  ///
  /// If T is an integer type (std::numeric_traits<T>::is_exact is
  /// true), then epsilon() returns 0.  Otherwise, if T is a
  /// floating-point type, it returns machine epsilon that T.
  static KOKKOS_FUNCTION mag_type epsilon();

  //@{
  /// \name Traits defined for backwards compatibility with
  /// Teuchos::ScalarTraits
  ///
  /// All of the typedefs, \c bool constants, and class methods in
  /// this section are defined in order that one may replace most uses
  /// of Teuchos::ScalarTraits with ArithTraits.  Users who do not
  /// have this backwards compatibility requirement should prefer
  /// equivalents in other sections.  Those class methods which have
  /// the same name and meaning in both Teuchos::ScalarTraits and this
  /// class, such as log() and pow(), are not in this section.

  //! Same as mag_type; the type of the absolute value (magnitude) of T.
  using magnitudeType = T;

  /// \brief The type with "half the precision" of T.
  ///
  /// This typedef only makes sense if T is a floating-point type.
  using halfPrecision = T;

  /// \brief The type with "twice the the precision" of T.
  ///
  /// This typedef only makes sense if T is a floating-point type.
  using doublePrecision = T;

  static constexpr bool isComplex    = false;
  static constexpr bool isOrdinal    = false;
  static constexpr bool isComparable = false;

  /// \brief True if this type T has floating-point parameters.
  ///
  /// This is true if and only if this specialization of ArithTraits
  /// has "machine-specific" parameters eps(), sfmin(), base(),
  /// prec(), t(), rnd(), emin(), rmin(), emax(), and rmax(), relating
  /// to floating-point types.
  static constexpr bool hasMachineParameters = false;

  //! Return relative machine precision.
  static KOKKOS_FUNCTION mag_type eps();

  //! Return safe minimum (sfmin), such that 1/sfmin does not overflow.
  static KOKKOS_FUNCTION mag_type sfmin();

  //! Return the base of the scalar type T.
  static KOKKOS_FUNCTION int base();

  //! Return <tt>eps*base</tt>.
  static KOKKOS_FUNCTION mag_type prec();

  //! Returns the number of (base) digits in the significand.
  static KOKKOS_FUNCTION int t();

  //! 1.0 when rounding occurs in addition, else 0.0.
  static KOKKOS_FUNCTION mag_type rnd();

  //! Returns the minimum exponent before (gradual) underflow.
  static KOKKOS_FUNCTION int emin();

  //! Returns the underflow threshold: <tt>base^(emin-1)</tt>
  static KOKKOS_FUNCTION mag_type rmin();

  //! Returns the largest exponent before overflow.
  static KOKKOS_FUNCTION int emax();

  //! Overflow theshold: <tt>(base^emax)*(1-eps)</tt>
  static KOKKOS_FUNCTION mag_type rmax();

  //! Same as abs(); return the magnitude of x.
  static KOKKOS_FUNCTION magnitudeType magnitude(const T& x);

  //! Same as conj(); return the complex conjugate of x.
  static KOKKOS_FUNCTION T conjugate(const T& x);

  /// \brief Whether x is (silent) NaN or Inf.
  ///
  /// This is the same as <tt>isNan(x) || isInf(x)</tt>.
  static KOKKOS_FUNCTION bool isnaninf(const T& x);

  /// \brief The string name of T.
  ///
  /// Note that this is not a device function.
  static std::string name();

  //! Same as sqrt(x); the square root of x.
  static KOKKOS_FUNCTION T squareroot(const T& x);
  //@}
};

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
template <>
class ArithTraits<Kokkos::Experimental::half_t> {
 public:
  using val_type = Kokkos::Experimental::half_t;
  using mag_type = val_type;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed      = true;
  static constexpr bool is_integer     = false;
  static constexpr bool is_exact       = false;
  static constexpr bool is_complex     = false;
  static constexpr bool has_infinity   = true;

  // Backwards compatibility with Teuchos::ScalarTraits.
  using magnitudeType   = mag_type;
  using halfPrecision   = Kokkos::Experimental::half_t;
  using doublePrecision = float;

  static std::string name() { return "half_t"; }

  static constexpr bool isComplex            = false;
  static constexpr bool isOrdinal            = false;
  static constexpr bool isComparable         = true;
  static constexpr bool hasMachineParameters = true;

#if defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_HIP)
  KOKKOSKERNELS_ARITHTRAITS_HALF_FP(KOKKOS_FUNCTION)
#else
  KOKKOSKERNELS_ARITHTRAITS_REAL_FP(KOKKOS_FUNCTION)
#endif
};
#endif  // #if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT

// Since Kokkos::Experimental::bhalf_t falls back to float, only define
// ArithTraits if bhalf_t is a backend specialization
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
template <>
class ArithTraits<Kokkos::Experimental::bhalf_t> {
 public:
  using val_type = Kokkos::Experimental::bhalf_t;
  using mag_type = val_type;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed      = true;
  static constexpr bool is_integer     = false;
  static constexpr bool is_exact       = false;
  static constexpr bool is_complex     = false;
  static constexpr bool has_infinity   = true;

  // Backwards compatibility with Teuchos::ScalarTraits.
  using magnitudeType  = mag_type;
  using bhalfPrecision = Kokkos::Experimental::bhalf_t;
  // There is no type that has twice the precision as bhalf_t.
  // The closest type would be float.
  using doublePrecision = void;

  static constexpr bool isComplex            = false;
  static constexpr bool isOrdinal            = false;
  static constexpr bool isComparable         = true;
  static constexpr bool hasMachineParameters = true;

  static std::string name() { return "bhalf_t"; }

#if defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_HIP)
  KOKKOSKERNELS_ARITHTRAITS_HALF_FP(KOKKOS_FUNCTION)
#else
  KOKKOSKERNELS_ARITHTRAITS_REAL_FP(KOKKOS_FUNCTION)
#endif
};
#endif  // #if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT

template <>
class ArithTraits<float> {
 public:
  using val_type = float;
  using mag_type = val_type;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed      = true;
  static constexpr bool is_integer     = false;
  static constexpr bool is_exact       = false;
  static constexpr bool is_complex     = false;
  static constexpr bool has_infinity   = true;

  // Backwards compatibility with Teuchos::ScalarTraits.
  using magnitudeType   = mag_type;
  using halfPrecision   = float;  // Should we switch to Kokkos::half_t
  using doublePrecision = double;

  static constexpr bool isComplex            = false;
  static constexpr bool isOrdinal            = false;
  static constexpr bool isComparable         = true;
  static constexpr bool hasMachineParameters = true;

  static std::string name() { return "float"; }

  KOKKOSKERNELS_ARITHTRAITS_REAL_FP(KOKKOS_FUNCTION)
};

template <>
class ArithTraits<double> {
 public:
  using val_type = double;
  using mag_type = val_type;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed      = true;
  static constexpr bool is_integer     = false;
  static constexpr bool is_exact       = false;
  static constexpr bool is_complex     = false;
  static constexpr bool has_infinity   = true;

  // Backwards compatibility with Teuchos::ScalarTraits.
  using magnitudeType = mag_type;
  using halfPrecision = float;
#if defined(__CUDA_ARCH__)
  using doublePrecision = double;  // CUDA doesn't support long double, unfortunately
#elif defined(__HIP_DEVICE_COMPILE__)
  using doublePrecision = double;  // HIP does not support long double unfortunately
#else
  using doublePrecision = long double;
#endif  // __CUDA_ARCH__
  static constexpr bool isComplex            = false;
  static constexpr bool isOrdinal            = false;
  static constexpr bool isComparable         = true;
  static constexpr bool hasMachineParameters = true;

  static std::string name() { return "double"; }

  KOKKOSKERNELS_ARITHTRAITS_REAL_FP(KOKKOS_FUNCTION)
};

// CUDA and HIP do not support long double in device functions,
// so none of the class methods in this specialization are marked
// as device functions.
template <>
class ArithTraits<long double> {
 public:
  using val_type = long double;
  using mag_type = long double;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed      = true;
  static constexpr bool is_integer     = false;
  static constexpr bool is_exact       = false;
  static constexpr bool is_complex     = false;
  static constexpr bool has_infinity   = true;

  // Backwards compatibility with Teuchos::ScalarTraits.
  using magnitudeType = mag_type;
  using halfPrecision = double;
  // It might be appropriate to use QD's qd_real here.
  // For now, long double is the most you get.
  using doublePrecision = val_type;

  static constexpr bool isComplex            = false;
  static constexpr bool isOrdinal            = false;
  static constexpr bool isComparable         = true;
  static constexpr bool hasMachineParameters = true;

  static std::string name() { return "long double"; }

  KOKKOSKERNELS_ARITHTRAITS_REAL_FP()
};  // long double specialization

#if defined(KOKKOS_ENABLE_LIBQUADMATH)
// CUDA does not support __float128 in device functions, so none of
// the class methods in this specialization are marked as device
// functions.
template <>
class ArithTraits<__float128> {
 public:
  using val_type = __float128;
  using mag_type = val_type;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed      = true;
  static constexpr bool is_integer     = false;
  static constexpr bool is_exact       = false;
  static constexpr bool is_complex     = false;
  static constexpr bool has_infinity   = true;

  // Backwards compatibility with Teuchos::ScalarTraits.
  using magnitudeType = mag_type;
  using halfPrecision = double;
  // Unfortunately, we can't rely on a standard __float256 type.
  using doublePrecision = __float128;

  static constexpr bool isComplex            = false;
  static constexpr bool isOrdinal            = false;
  static constexpr bool isComparable         = true;
  static constexpr bool hasMachineParameters = true;

  static std::string name() { return "__float128"; }

  KOKKOSKERNELS_ARITHTRAITS_REAL_FP()
};      // __float128 specialization
#endif  // KOKKOS_ENABLE_LIBQUADMATH

template <>
class ArithTraits< ::Kokkos::complex<float> > {
 public:
  using val_type = ::Kokkos::complex<float>;
  using mag_type = float;

  static std::string name() { return "Kokkos::complex<float>"; }

  KOKKOSKERNELS_ARITHTRAITS_CMPLX_FP(KOKKOS_FUNCTION)
};

template <>
class ArithTraits< ::Kokkos::complex<double> > {
 public:
  using val_type = ::Kokkos::complex<double>;
  using mag_type = double;

  static std::string name() { return "Kokkos::complex<double>"; }

  KOKKOSKERNELS_ARITHTRAITS_CMPLX_FP(KOKKOS_FUNCTION)
};

/// \brief Partial specialization for std::complex<RealFloatType>.
///
/// The C++ Standard Library (with C++03 at least) only allows
/// std::complex<RealFloatType> for RealFloatType = float, double, or
/// long double.
template <class RealFloatType>
class ArithTraits<std::complex<RealFloatType> > {
 public:
  //! Kokkos internally replaces std::complex with Kokkos::complex.
  using val_type = ::Kokkos::complex<RealFloatType>;
  using mag_type = RealFloatType;

  static constexpr bool is_specialized = true;
  static constexpr bool is_signed      = true;
  static constexpr bool is_integer     = false;
  static constexpr bool is_exact       = false;
  static constexpr bool is_complex     = true;

  static constexpr bool has_infinity = true;
  static std::complex<RealFloatType> infinity() {
    return std::complex<RealFloatType>(ArithTraits<mag_type>::infinity(), ArithTraits<mag_type>::infinity());
  }

#ifdef KOKKOS_ENABLE_SYCL
  template <typename Dummy = RealFloatType>
  static bool isInf(const std::complex<Dummy>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isinf;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL)
    using sycl::isinf;
#endif
    return isinf(real(x)) || isinf(imag(x));
  }
  template <>
  static bool isInf<long double>(const std::complex<long double>& x) {
    Kokkos::abort("isInf not available for std::complex<long double>!\n");
    return true;
  }
#else
  static bool isInf(const std::complex<RealFloatType>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isinf;
#endif
    return isinf(real(x)) || isinf(imag(x));
  }
#endif
#ifdef KOKKOS_ENABLE_SYCL
  template <typename Dummy = RealFloatType>
  static bool isNan(const std::complex<Dummy>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isnan;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL)
    using sycl::isnan;
#endif
    return isnan(real(x)) || isnan(imag(x));
  }
  template <>
  static bool isNan<long double>(const std::complex<long double>& x) {
    Kokkos::abort("isNan not available for std::complex<long double>!\n");
    return true;
  }
#else
  static bool isNan(const std::complex<RealFloatType>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isnan;
#endif
    return isnan(real(x)) || isnan(imag(x));
  }
#endif
  static mag_type abs(const std::complex<RealFloatType>& x) { return std::abs(x); }
  static std::complex<RealFloatType> zero() {
    return std::complex<RealFloatType>(ArithTraits<mag_type>::zero(), ArithTraits<mag_type>::zero());
  }
  static std::complex<RealFloatType> one() {
    return std::complex<RealFloatType>(ArithTraits<mag_type>::one(), ArithTraits<mag_type>::zero());
  }
  static std::complex<RealFloatType> min() {
    return std::complex<RealFloatType>(ArithTraits<mag_type>::min(), ArithTraits<mag_type>::zero());
  }
  static std::complex<RealFloatType> max() {
    return std::complex<RealFloatType>(ArithTraits<mag_type>::max(), ArithTraits<mag_type>::zero());
  }
  static mag_type real(const std::complex<RealFloatType>& x) { return std::real(x); }
  static mag_type imag(const std::complex<RealFloatType>& x) { return std::imag(x); }
  static std::complex<RealFloatType> conj(const std::complex<RealFloatType>& x) { return std::conj(x); }
  static std::complex<RealFloatType> pow(const std::complex<RealFloatType>& x, const std::complex<RealFloatType>& y) {
    // Fix for some weird gcc 4.2.1 inaccuracy.
    if (y == one()) {
      return x;
    } else if (y == one() + one()) {
      return x * x;
    } else {
      return std::pow(x, y);
    }
  }
  static std::complex<RealFloatType> pow(const std::complex<RealFloatType>& x, const RealFloatType& y) {
    // Fix for some weird gcc 4.2.1 inaccuracy.
    if (y == ArithTraits<RealFloatType>::one()) {
      return x;
    } else if (y == ArithTraits<RealFloatType>::one() + ArithTraits<RealFloatType>::one()) {
      return x * x;
    } else {
      return std::pow(x, y);
    }
  }
  static std::complex<RealFloatType> sqrt(const std::complex<RealFloatType>& x) { return std::sqrt(x); }
  static std::complex<RealFloatType> cbrt(const std::complex<RealFloatType>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::cbrt(x);
#else
    return ::cbrt(x);
#endif
  }
  static std::complex<RealFloatType> exp(const std::complex<RealFloatType>& x) { return std::exp(x); }
  static std::complex<RealFloatType> log(const std::complex<RealFloatType>& x) { return std::log(x); }
  static std::complex<RealFloatType> log10(const std::complex<RealFloatType>& x) { return std::log10(x); }
  static std::complex<RealFloatType> sin(const std::complex<RealFloatType>& x) { return std::sin(x); }
  static std::complex<RealFloatType> cos(const std::complex<RealFloatType>& x) { return std::cos(x); }
  static std::complex<RealFloatType> tan(const std::complex<RealFloatType>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::tan(x);
#else
    return std::tan(x);
#endif
  }
  static std::complex<RealFloatType> sinh(const std::complex<RealFloatType>& x) { return std::sinh(x); }
  static std::complex<RealFloatType> cosh(const std::complex<RealFloatType>& x) { return std::cosh(x); }
  static std::complex<RealFloatType> tanh(const std::complex<RealFloatType>& x) { return std::tanh(x); }
  static std::complex<RealFloatType> asin(const std::complex<RealFloatType>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::asin(x);
#else
    return ::asin(x);
#endif
  }
  static std::complex<RealFloatType> acos(const std::complex<RealFloatType>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::acos(x);
#else
    return ::acos(x);
#endif
  }
  static std::complex<RealFloatType> atan(const std::complex<RealFloatType>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    using sycl::atan;
#else
    using std::atan;
#endif
    return atan(x);
  }
  static std::complex<RealFloatType> nan() {
    const mag_type mag_nan = ArithTraits<mag_type>::nan();
    return std::complex<RealFloatType>(mag_nan, mag_nan);
  }
  static mag_type epsilon() { return ArithTraits<mag_type>::epsilon(); }

  // Backwards compatibility with Teuchos::ScalarTraits.
  using magnitudeType   = mag_type;
  using halfPrecision   = std::complex<typename ArithTraits<mag_type>::halfPrecision>;
  using doublePrecision = std::complex<typename ArithTraits<mag_type>::doublePrecision>;

  static constexpr bool isComplex            = true;
  static constexpr bool isOrdinal            = false;
  static constexpr bool isComparable         = false;
  static constexpr bool hasMachineParameters = true;
  static bool isnaninf(const std::complex<RealFloatType>& x) { return isNan(x) || isInf(x); }
  static mag_type magnitude(const std::complex<RealFloatType>& x) { return abs(x); }
  static std::complex<RealFloatType> conjugate(const std::complex<RealFloatType>& x) { return conj(x); }
  static std::string name() { return std::string("std::complex<") + ArithTraits<mag_type>::name() + ">"; }
  static std::complex<RealFloatType> squareroot(const std::complex<RealFloatType>& x) { return sqrt(x); }
  static mag_type eps() { return epsilon(); }
  static mag_type sfmin() { return ArithTraits<mag_type>::sfmin(); }
  static int base() { return ArithTraits<mag_type>::base(); }
  static mag_type prec() { return ArithTraits<mag_type>::prec(); }
  static int t() { return ArithTraits<mag_type>::t(); }
  static mag_type rnd() { return ArithTraits<mag_type>::one(); }
  static int emin() { return ArithTraits<mag_type>::emin(); }
  static mag_type rmin() { return ArithTraits<mag_type>::rmin(); }
  static int emax() { return ArithTraits<mag_type>::emax(); }
  static mag_type rmax() { return ArithTraits<mag_type>::rmax(); }
};

template <>
class ArithTraits<char> {
 public:
  using val_type = char;
  using mag_type = val_type;

  // The C(++) standard does not require that char be signed.  In
  // fact, signed char, unsigned char, and char are distinct types.
  // We can use std::numeric_limits here because it's a const bool,
  // not a class method.
  static constexpr bool is_signed = std::numeric_limits<val_type>::is_signed;

  static std::string name() { return "char"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<signed char> {
 public:
  using val_type = signed char;
  using mag_type = val_type;

  static constexpr bool is_signed = true;

  static std::string name() { return "signed char"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<unsigned char> {
 public:
  using val_type = unsigned char;
  using mag_type = val_type;

  static constexpr bool is_signed = false;

  static std::string name() { return "unsigned char"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<short> {
 public:
  using val_type = short;
  using mag_type = val_type;

  static constexpr bool is_signed = true;

  static std::string name() { return "short"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<unsigned short> {
 public:
  using val_type = unsigned short;
  using mag_type = val_type;

  static constexpr bool is_signed = false;

  static std::string name() { return "unsigned short"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<int> {
 public:
  using val_type = int;
  using mag_type = val_type;

  static constexpr bool is_signed = true;

  static std::string name() { return "int"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<unsigned int> {
 public:
  using val_type = unsigned int;
  using mag_type = val_type;

  static constexpr bool is_signed = false;

  static std::string name() { return "unsigned int"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<long> {
 public:
  using val_type = long;
  using mag_type = val_type;

  static constexpr bool is_signed = true;

  static std::string name() { return "long"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<unsigned long> {
 public:
  using val_type = unsigned long;
  using mag_type = val_type;

  static constexpr bool is_signed = false;

  static std::string name() { return "unsigned long"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<long long> {
 public:
  using val_type = long long;
  using mag_type = val_type;

  static constexpr bool is_signed = true;

  static std::string name() { return "long long"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

template <>
class ArithTraits<unsigned long long> {
 public:
  using val_type = unsigned long long;
  using mag_type = val_type;

  static constexpr bool is_signed = false;

  static std::string name() { return "unsigned long long"; }

  KOKKOSKERNELS_ARITHTRAITS_INTEGRAL()
};

// dd_real and qd_real are floating-point types provided by the QD
// library of David Bailey (LBNL):
//
// http://crd-legacy.lbl.gov/~dhbailey/mpdist/
//
// dd_real uses two doubles (128 bits), and qd_real uses four doubles
// (256 bits).
//
// Kokkos does <i>not</i> currently support these types in device
// functions.  It should be possible to use Kokkos' support for
// aggregate types to implement device function support for dd_real
// and qd_real, but we have not done this yet (as of 09 Jan 2015).
// Hence, the class methods of the ArithTraits specializations for
// dd_real and qd_real are not marked as device functions.
#ifdef HAVE_KOKKOS_QD
// LBV: I would like to deprecate this strange optional
// dependency on the lbnl package, is there anyone actully
// using this? It certainly is never tested by CI or nightly
// so probably does not work...
template <>
struct [[deprecated]] ArithTraits<dd_real> {
  typedef dd_real val_type;
  typedef dd_real mag_type;

  static const bool is_specialized = true;
  static const bool is_signed      = true;
  static const bool is_integer     = false;
  static const bool is_exact       = false;
  static const bool is_complex     = false;

  static inline bool isInf(const val_type& x) { return isinf(x); }
  static inline bool isNan(const val_type& x) { return isnan(x); }
  static inline mag_type abs(const val_type& x) { return ::abs(x); }
  static inline val_type zero() { return val_type(0.0); }
  static inline val_type one() { return val_type(1.0); }
  static inline val_type min() { return std::numeric_limits<val_type>::min(); }
  static inline val_type max() { return std::numeric_limits<val_type>::max(); }
  static inline mag_type real(const val_type& x) { return x; }
  static inline mag_type imag(const val_type&) { return zero(); }
  static inline val_type conj(const val_type& x) { return x; }
  static inline val_type pow(const val_type& x, const val_type& y) { return ::pow(x, y); }
  static inline val_type sqrt(const val_type& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::sqrt(x);
#else
    return ::sqrt(x);
#endif
  }
  static inline val_type cbrt(const val_type& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::cbrt(x);
#else
    return ::cbrt(x);
#endif
  }
  static inline val_type exp(const val_type& x) { return ::exp(x); }
  static inline val_type log(const val_type& x) {
    // dd_real puts its transcendental functions in the global namespace.
    return ::log(x);
  }
  static inline val_type log10(const val_type& x) { return ::log10(x); }
  static KOKKOS_FUNCTION val_type sin(const val_type x) { return ::sin(x); }
  static KOKKOS_FUNCTION val_type cos(const val_type x) { return ::cos(x); }
  static KOKKOS_FUNCTION val_type tan(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::tan(x);
#else
    return std::tan(x);
#endif
  }
  static KOKKOS_FUNCTION val_type sinh(const val_type x) { return ::sinh(x); }
  static KOKKOS_FUNCTION val_type cosh(const val_type x) { return ::cosh(x); }
  static KOKKOS_FUNCTION val_type tanh(const val_type x) { return ::tanh(x); }
  static KOKKOS_FUNCTION val_type asin(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::asin(x);
#else
    return ::asin(x);
#endif
  }
  static KOKKOS_FUNCTION val_type acos(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::acos(x);
#else
    return ::acos(x);
#endif
  }
  static KOKKOS_FUNCTION val_type atan(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::atan(x);
#else
    return ::atan(x);
#endif
  }
  static inline val_type nan() { return val_type::_nan; }
  static val_type epsilon() { return std::numeric_limits<val_type>::epsilon(); }

  typedef dd_real magnitudeType;
  typedef double halfPrecision;
  typedef qd_real doublePrecision;

  static const bool isComplex            = false;
  static const bool isOrdinal            = false;
  static const bool isComparable         = true;
  static const bool hasMachineParameters = true;

  static mag_type eps() { return epsilon(); }
  static mag_type sfmin() { return min(); }
  static int base() { return std::numeric_limits<val_type>::radix; }
  static mag_type prec() { return eps() * base(); }
  static int t() { return std::numeric_limits<val_type>::digits; }
  static mag_type rnd() { return std::numeric_limits<val_type>::round_style == std::round_to_nearest ? one() : zero(); }
  static int emin() { return std::numeric_limits<val_type>::min_exponent; }
  static mag_type rmin() { return std::numeric_limits<val_type>::min(); }
  static int emax() { return std::numeric_limits<val_type>::max_exponent; }
  static mag_type rmax() { return std::numeric_limits<val_type>::max(); }
  static mag_type magnitude(const val_type& x) { return ::abs(x); }
  static val_type conjugate(const val_type& x) { return conj(x); }
  static bool isnaninf(const val_type& x) { return isNan(x) || isInf(x); }
  static std::string name() { return "dd_real"; }
  static val_type squareroot(const val_type& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::sqrt(x);
#else
    return ::sqrt(x);
#endif
  }
};

template <>
struct [[deprecated]] ArithTraits<qd_real> {
  typedef qd_real val_type;
  typedef qd_real mag_type;

  static const bool is_specialized = true;
  static const bool is_signed      = true;
  static const bool is_integer     = false;
  static const bool is_exact       = false;
  static const bool is_complex     = false;

  static inline bool isInf(const val_type& x) { return isinf(x); }
  static inline bool isNan(const val_type& x) { return isnan(x); }
  static inline mag_type abs(const val_type& x) { return ::abs(x); }
  static inline val_type zero() { return val_type(0.0); }
  static inline val_type one() { return val_type(1.0); }
  static inline val_type min() { return std::numeric_limits<val_type>::min(); }
  static inline val_type max() { return std::numeric_limits<val_type>::max(); }
  static inline mag_type real(const val_type& x) { return x; }
  static inline mag_type imag(const val_type&) { return zero(); }
  static inline val_type conj(const val_type& x) { return x; }
  static inline val_type pow(const val_type& x, const val_type& y) { return ::pow(x, y); }
  static inline val_type sqrt(const val_type& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::sqrt(x);
#else
    return ::sqrt(x);
#endif
  }
  static inline val_type cbrt(const val_type& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::cbrt(x);
#else
    return ::cbrt(x);
#endif
  }
  static inline val_type exp(const val_type& x) { return ::exp(x); }
  static inline val_type log(const val_type& x) {
    // val_type puts its transcendental functions in the global namespace.
    return ::log(x);
  }
  static inline val_type log10(const val_type& x) { return ::log10(x); }
  static KOKKOS_FUNCTION val_type sin(const val_type x) { return ::sin(x); }
  static KOKKOS_FUNCTION val_type cos(const val_type x) { return ::cos(x); }
  static KOKKOS_FUNCTION val_type tan(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::tan(x);
#else
    return std::tan(x);
#endif
  }
  static KOKKOS_FUNCTION val_type sinh(const val_type x) { return ::sinh(x); }
  static KOKKOS_FUNCTION val_type cosh(const val_type x) { return ::cosh(x); }
  static KOKKOS_FUNCTION val_type tanh(const val_type x) { return ::tanh(x); }
  static KOKKOS_FUNCTION val_type asin(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::asin(x);
#else
    return ::asin(x);
#endif
  }
  static KOKKOS_FUNCTION val_type acos(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::acos(x);
#else
    return ::acos(x);
#endif
  }
  static KOKKOS_FUNCTION val_type atan(const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::atan(x);
#else
    return ::atan(x);
#endif
  }
  static inline val_type nan() { return val_type::_nan; }
  static inline val_type epsilon() { return std::numeric_limits<val_type>::epsilon(); }

  typedef qd_real magnitudeType;
  typedef dd_real halfPrecision;
  // The QD library does not have an "oct-double real" class.  One
  // could use an arbitrary-precision library like MPFR or ARPREC,
  // with the precision set appropriately, to get an
  // extended-precision type for qd_real.
  typedef qd_real doublePrecision;

  static const bool isComplex            = false;
  static const bool isOrdinal            = false;
  static const bool isComparable         = true;
  static const bool hasMachineParameters = true;

  static mag_type eps() { return epsilon(); }
  static mag_type sfmin() { return min(); }
  static int base() { return std::numeric_limits<val_type>::radix; }
  static mag_type prec() { return eps() * base(); }
  static int t() { return std::numeric_limits<val_type>::digits; }
  static mag_type rnd() { return std::numeric_limits<val_type>::round_style == std::round_to_nearest ? one() : zero(); }
  static int emin() { return std::numeric_limits<val_type>::min_exponent; }
  static mag_type rmin() { return std::numeric_limits<val_type>::min(); }
  static int emax() { return std::numeric_limits<val_type>::max_exponent; }
  static mag_type rmax() { return std::numeric_limits<val_type>::max(); }
  static mag_type magnitude(const val_type& x) { return ::abs(x); }
  static val_type conjugate(const val_type& x) { return conj(x); }
  static bool isnaninf(const val_type& x) { return isNan(x) || isInf(x); }
  static std::string name() { return "qd_real"; }
  static val_type squareroot(const val_type& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    return sycl::sqrt(x);
#else
    return ::sqrt(x);
#endif
  }
};
#endif  // HAVE_KOKKOS_QD

namespace Details {
template <typename T>
using ArithTraits [[deprecated("Use Kokkos::ArithTraits instead")]] = ::Kokkos::ArithTraits<T>;

}  // namespace Details
}  // namespace Kokkos

#endif  // KOKKOS_ARITHTRAITS_HPP
