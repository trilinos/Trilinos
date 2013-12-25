/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_ARITHTRAITS_HPP
#define KOKKOS_ARITHTRAITS_HPP

/// \file Kokkos_ArithTraits.hpp
/// \brief Declaration and definition of Kokkos::Details::ArithTraits
///
/// \warning This interface is not ready for exposure to users yet.
///   Users beware!
///
/// This file works with CUDA 5.5 (with gcc, see next), gcc 4.2.1
/// (stock Mac XCode), and Clang 3.2, all on Mac.
///
/// On my Mac, for all of the above compilers, long and int64_t are
/// different, even though sizeof(long) == 8.  This manifests by an
/// int64_t specialization not sufficing for long.  This could be
/// because long and long long are two different types (not aliases of
/// one another), even though they might have the same size.
/// Interestingly, though, int64_t and long long do appear to be the
/// same type on my system.

#include <cfloat>
#include <climits>
#include <cmath>
#include <stdint.h> // CUDA 5.5 doesn't recognize <cstdint>.
#include <complex> // std::complex

//
// mfh 24 Dec 2013: Temporary measure for testing; will go away.
//
#ifndef KOKKOS_DEVICE_FUNCTION
#  ifdef __CUDA_ARCH__
#    define KOKKOS_DEVICE_FUNCTION inline __host__ __device__
#  else
#    define KOKKOS_DEVICE_FUNCTION
#  endif // __CUDA_ARCH__
#endif // KOKKOS_DEVICE_FUNCTION

namespace { // anonymous

/// \fn intPowImpl
/// \tparam IntType A built-in integer type.
/// \brief Implementation of intPowSigned and intPowUnsigned.
///
/// \pre x != 0
/// \pre y > 0
///
/// Use intPowSigned or intPowUnsigned for general y.
template<class IntType>
KOKKOS_DEVICE_FUNCTION IntType
intPowImpl (const IntType x, const IntType y)
{
  // Recursion (unrolled into while loop): pow(x, 2y) = (x^y)^2
  IntType prod = x;
  IntType y_cur = 1;
  // If y == 1, then prod stays x.
  while (y_cur < y) {
    prod = prod * prod;
    y_cur = y_cur << 1;
  }
  // abs(y - y_cur) < floor(log2(y)), so it won't hurt asymptotic run
  // time to finish the remainder in a linear iteration.
  if (y > y_cur) {
    const IntType left = y - y_cur;
    for (IntType k = 0; k < left; ++k) {
      prod = prod * x;
    }
  }
  else if (y < y_cur) {
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



/// \fn intPowSigned
/// \tparam IntType A built-in signed integer type.
/// \brief Compute x raised to the power y.
///
/// If the arguments are invalid (e.g., if x and y are both zero), the
/// result of this function is undefined.  However, this function will
/// not throw an exception in that case.
template<class IntType>
KOKKOS_DEVICE_FUNCTION IntType
intPowSigned (const IntType x, const IntType y)
{
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
    }
    else if (x == -1) {
      return (y % 2 == 0) ? 1 : -1;
    }
    else {
      return 0; // round the fraction to zero
    }
  } else {
    return intPowImpl<IntType> (x, y);
  }
}

/// \fn intPowUnsigned
/// \tparam IntType A built-in unsigned integer type.
/// \brief Compute x raised to the power y.
///
/// If the arguments are invalid (e.g., if x and y are both zero), the
/// result of this function is undefined.  However, this function will
/// not throw an exception in that case.
template<class IntType>
KOKKOS_DEVICE_FUNCTION IntType
intPowUnsigned (const IntType x, const IntType y)
{
  // It's not entirely clear what to return if x and y are both zero.
  // In the case of floating-point numbers, 0^0 is NaN.  Here, though,
  // I think it's safe to return 0.
  if (x == 0) {
    return 0;
  } else if (y == 0) {
    return 1;
  } else {
    return intPowImpl<IntType> (x, y);
  }
}

// // http://www.azillionmonkeys.com/qed/sqroot.html#implementations
//
// // This is not technically correct, because ANSI C(++) aliasing
// // rules forbid things like the assignment to tempf below.
// int32_t isqrt (const int32_t r) {
//   float tempf, x, y, rr;
//   int32_t is;
//
//   rr = (int32_t) r;
//   y = rr * 0.5;
//   *(uint32_t*) &tempf = (0xbe6f0000 - *(uint32_t*) &rr) >> 1;
//   x = tempf;
//   x = (1.5*x) - (x*x)*(x*y);
//   if (r > 101123) {
//     x = (1.5*x) - (x*x)*(x*y);
//   }
//   is = (int32_t) (x*rr + 0.5);
//   return is + ((int32_t) (r - is*is)) >> 31;
// }

} // namespace (anonymous)

namespace Kokkos {
namespace Details {

/// \class ArithTraits
/// \brief Traits class for arithmetic on type T.
/// \tparam T "Scalar" type of interest
///
/// \warning This interface is not ready for exposure to users yet.
///   Users beware!
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
/// using the KOKKOS_DEVICE_FUNCTION macro.  All class methods must be
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
/// \section Kokkos_ArithTraits_unsupp Unsupported types
///
/// CUDA does not support long double or std::complex<T> in device
/// functions.  ArithTraits does have specializations for these types,
/// but the class methods therein are not marked as device functions.
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
template<class T>
class ArithTraits {
public:
  //! The type T itself.
  typedef T val_type;
  //! The type of the magnitude (absolute value) of T.
  typedef T mag_type;

  //! Whether ArithTraits has a specialization for T.
  static const bool is_specialized = false;
  //! Whether T is a signed type (has negative values).
  static const bool is_signed = false;
  //! Whether T is an integer type.
  static const bool is_integer = false;
  //! Whether T "uses exact representations" (or is approximate, i.e., may commit rounding error).
  static const bool is_exact = false;
  //! Whether T is a complex-valued type.
  static const bool is_complex = false;

  /// \brief Whether x is Inf.
  ///
  /// This can only be true for floating-point types T that support
  /// Inf.  If T is a complex type, we say that a T instance x is Inf
  /// if and only if <tt>isinf(real(x)) || isinf(imag(x))</tt>.
  ///
  /// Unfortunately we can't call this "isinf" (the equivalent C99
  /// function), because CUDA appears to implement that function using
  /// a macro, rather than using a function (as C++11 requires).
  static KOKKOS_DEVICE_FUNCTION bool isInf (const T& x);

  /// \brief Whether x is NaN (not a number).
  ///
  /// This can only be true for floating-point types T that support
  /// NaN.  If T is a complex type, we say that a T instance x is NaN
  /// if and only if <tt>isNan(real(x)) || isNan(imag(x))</tt>.
  ///
  /// Unfortunately we can't call this "isnan" (the equivalent C99
  /// function), because CUDA appears to implement that function using
  /// a macro, rather than using a function (as C++11 requires).
  static KOKKOS_DEVICE_FUNCTION bool isNan (const T& x);

  //! The absolute value (magnitude) of x.
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const T& x);

  //! The zero value of T; the arithmetic identity.
  static KOKKOS_DEVICE_FUNCTION T zero ();

  //! The one value of T; the multiplicative identity.
  static KOKKOS_DEVICE_FUNCTION T one ();

  /// \brief The minimum possible value of T.
  ///
  /// If T is a real floating-point type, then this is the minimum
  /// <i>positive</i> value, as with std::numeric_limits<T>::min().
  static KOKKOS_DEVICE_FUNCTION T min ();

  //! The maximum possible value of T.
  static KOKKOS_DEVICE_FUNCTION T max ();

  /// \brief The real part of x.
  ///
  /// If \c is_complex is false, then this just returns x.
  static KOKKOS_DEVICE_FUNCTION mag_type real (const T& x);

  /// \brief The imaginary part of x.
  ///
  /// If \c is_complex is false, then this just returns zero().
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const T&);

  /// \brief The complex conjugate of x.
  ///
  /// If \c is_complex is false, then this just returns x.
  static KOKKOS_DEVICE_FUNCTION T conj (const T&);

  //! x raised to the power y.
  static KOKKOS_DEVICE_FUNCTION T pow (const T& x, const T& y);

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
  static KOKKOS_DEVICE_FUNCTION T sqrt (const T& x);

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
  static KOKKOS_DEVICE_FUNCTION T log (const T& x);

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
  static KOKKOS_DEVICE_FUNCTION T log10 (const T& x);

  /// \brief Return a silent NaN, if appropriate for T.
  ///
  /// If T does <i>not</i> implement a silent NaN, the return value is
  /// undefined, but calling this method is still allowed.
  static KOKKOS_DEVICE_FUNCTION T nan ();

  /// \brief Machine epsilon.
  ///
  /// If T is an integer type (std::numeric_traits<T>::is_exact is
  /// true), then epsilon() returns 0.  Otherwise, if T is a
  /// floating-point type, it returns machine epsilon that T.
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon ();

  //@{
  /// \name Traits defined for backwards compatibility with Teuchos::ScalarTraits
  ///
  /// All of the typedefs, \c bool constants, and class methods in
  /// this section are defined in order that one may replace most uses
  /// of Teuchos::ScalarTraits with ArithTraits.  Users who do not
  /// have this backwards compatibility requirement should prefer
  /// equivalents in other sections.  Those class methods which have
  /// the same name and meaning in both Teuchos::ScalarTraits and this
  /// class, such as log() and pow(), are not in this section.

  //! Same as mag_type; the type of the absolute value (magnitude) of T.
  typedef T magnitudeType;

  /// \brief The type with "half the precision" of T.
  ///
  /// This typedef only makes sense if T is a floating-point type.
  typedef T halfPrecision;

  /// \brief The type with "twice the the precision" of T.
  ///
  /// This typedef only makes sense if T is a floating-point type.
  typedef T doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = false;

  /// \brief True if this type T has floating-point parameters.
  ///
  /// This is true if and only if this specialization of ArithTraits
  /// has "machine-specific" parameters eps(), sfmin(), base(),
  /// prec(), t(), rnd(), emin(), rmin(), emax(), and rmax(), relating
  /// to floating-point types.
  static const bool hasMachineParameters = false;

  //! Return relative machine precision.
  static KOKKOS_DEVICE_FUNCTION mag_type eps ();

  //! Return safe minimum (sfmin), such that 1/sfmin does not overflow.
  static KOKKOS_DEVICE_FUNCTION mag_type sfmin ();

  //! Return the base of the scalar type T.
  static KOKKOS_DEVICE_FUNCTION int base ();

  //! Return <tt>eps*base</tt>.
  static KOKKOS_DEVICE_FUNCTION mag_type prec ();

  //! Returns the number of (base) digits in the significand.
  static KOKKOS_DEVICE_FUNCTION int t ();

  //! 1.0 when rounding occurs in addition, else 0.0.
  static KOKKOS_DEVICE_FUNCTION mag_type rnd ();

  //! Returns the minimum exponent before (gradual) underflow.
  static KOKKOS_DEVICE_FUNCTION int emin ();

  //! Returns the underflow threshold: <tt>base^(emin-1)</tt>
  static KOKKOS_DEVICE_FUNCTION mag_type rmin ();

  //! Returns the largest exponent before overflow.
  static KOKKOS_DEVICE_FUNCTION int emax ();

  //! Overflow theshold: <tt>(base^emax)*(1-eps)</tt>
  static KOKKOS_DEVICE_FUNCTION mag_type rmax ();

  //! Same as abs(); return the magnitude of x.
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const T& x);

  //! Same as conj(); return the complex conjugate of x.
  static KOKKOS_DEVICE_FUNCTION T conjugate (const T& x);

  /// \brief Whether x is (silent) NaN or Inf.
  ///
  /// This is the same as <tt>isNan(x) || isInf(x)</tt>.
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const T& x);

  /// \brief The string name of T.
  ///
  /// Note that this is not a device function.
  static std::string name ();

  //! Same as sqrt(x); the square root of x.
  static KOKKOS_DEVICE_FUNCTION T squareroot (const T& x);
  //@}
};


template<>
class ArithTraits<float> {
public:
  typedef float val_type;
  typedef float mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const float x) {
#ifdef __CUDACC__
    return isinf (x);
#else
    return std::isinf (x);
#endif // __CUDACC__
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const float x) {
#ifdef __CUDACC__
    return isnan (x);
#else
    return std::isnan (x);
#endif // __CUDACC__
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const float x) {
    return ::fabs (x);
  }
  static KOKKOS_DEVICE_FUNCTION float zero () {
    return 0.0;
  }
  static KOKKOS_DEVICE_FUNCTION float one () {
    return 1.0;
  }
  static KOKKOS_DEVICE_FUNCTION float min () {
    return FLT_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION float max () {
    return FLT_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const float x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const float) {
    return 0.0;
  }
  static KOKKOS_DEVICE_FUNCTION float conj (const float x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION float pow (const float x, const float y) {
    return ::pow (x, y);
  }
  static KOKKOS_DEVICE_FUNCTION float sqrt (const float x) {
    return ::sqrt (x);
  }
  static KOKKOS_DEVICE_FUNCTION float log (const float x) {
    return ::log (x);
  }
  static KOKKOS_DEVICE_FUNCTION float log10 (const float x) {
    return ::log10 (x);
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return FLT_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef float magnitudeType;
  typedef float halfPrecision;
  typedef double doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const float x) {
    return isNan (x) || isInf (x);
  }
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const float x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION float conjugate (const float x) {
    return conj (x);
  }
  //static KOKKOS_DEVICE_FUNCTION bool isnaninf (const float);
  static std::string name () {
    return "float";
  }
  static KOKKOS_DEVICE_FUNCTION float squareroot (const float x) {
    return sqrt (x);
  }
  static KOKKOS_DEVICE_FUNCTION float nan () {
#ifdef __CUDA_ARCH__
    return nan ();
#else
    // http://pubs.opengroup.org/onlinepubs/009696899/functions/nan.html
    return strtof ("NAN()", (char**) NULL);
#endif // __CUDA_ARCH__
  }
  static KOKKOS_DEVICE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_DEVICE_FUNCTION mag_type sfmin () {
    return FLT_MIN; // ???
  }
  static KOKKOS_DEVICE_FUNCTION int base () {
    return FLT_RADIX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type prec () {
    return eps () * static_cast<mag_type> (base ());
  }
  static KOKKOS_DEVICE_FUNCTION int t () {
    return FLT_MANT_DIG;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type rnd () {
    return 1.0;
  }
  static KOKKOS_DEVICE_FUNCTION int emin () {
    return FLT_MIN_EXP;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type rmin () {
    return FLT_MIN; // ??? // should be base^(emin-1)
  }
  static KOKKOS_DEVICE_FUNCTION int emax () {
    return FLT_MAX_EXP;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type rmax () {
    return FLT_MAX; // ??? // should be (base^emax)*(1-eps)
  }
};



template<class RealFloatType>
class ArithTraits<std::complex<RealFloatType> > {
public:
  typedef std::complex<RealFloatType> val_type;
  typedef RealFloatType mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = true;

  static bool isInf (const val_type& x) {
#ifdef __CUDACC__
    return isinf (real (x)) || isinf (imag (x));
#else
    return std::isinf (real (x)) || std::isinf (imag (x));
#endif // __CUDACC__
  }
  static bool isNan (const val_type& x) {
#ifdef __CUDACC__
    return isnan (real (x)) || isnan (imag (x));
#else
    return std::isnan (real (x)) || std::isnan (imag (x));
#endif // __CUDACC__
  }
  static mag_type abs (const val_type& x) {
    return std::abs (x);
  }
  static val_type zero () {
    return val_type (ArithTraits<mag_type>::zero (), ArithTraits<mag_type>::zero ());
  }
  static val_type one () {
    return val_type (ArithTraits<mag_type>::one (), ArithTraits<mag_type>::zero ());
  }
  static val_type min () {
    return val_type (ArithTraits<mag_type>::min (), ArithTraits<mag_type>::zero ());
  }
  static val_type max () {
    return val_type (ArithTraits<mag_type>::max (), ArithTraits<mag_type>::zero ());
  }
  static mag_type real (const val_type& x) {
    return std::real (x);
  }
  static mag_type imag (const val_type& x) {
    return std::imag (x);
  }
  static val_type conj (const val_type& x) {
    return std::conj (x);
  }
  static val_type
  pow (const val_type& x, const val_type& y) {
    // Fix for some weird gcc 4.2.1 inaccuracy.
    if (y == one ()) {
      return x;
    } else if (y == one () + one ()) {
      return x * x;
    } else {
      return std::pow (x, y);
    }
  }
  static val_type sqrt (const val_type& x) {
    return std::sqrt (x);
  }
  static val_type log (const val_type& x) {
    return std::log (x);
  }
  static val_type log10 (const val_type& x) {
    return std::log10 (x);
  }
  static val_type nan () {
    const mag_type mag_nan = ArithTraits<mag_type>::nan ();
    return val_type (mag_nan, mag_nan);
  }
  static mag_type epsilon () {
    return ArithTraits<mag_type>::epsilon ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef std::complex<typename ArithTraits<mag_type>::halfPrecision> halfPrecision;
  typedef std::complex<typename ArithTraits<mag_type>::doublePrecision> doublePrecision;

  static const bool isComplex = true;
  static const bool isOrdinal = false;
  static const bool isComparable = false;
  static const bool hasMachineParameters = true;
  static bool isnaninf (const val_type& x) {
    return isNan (x) || isInf (x);
  }
  static mag_type magnitude (const val_type& x) {
    return abs (x);
  }
  static val_type conjugate (const val_type& x) {
    return conj (x);
  }
  //static KOKKOS_DEVICE_FUNCTION bool isnaninf (const float);
  static std::string name () {
    return std::string ("std::complex<") + ArithTraits<mag_type>::name () + ">";
  }
  static val_type squareroot (const val_type& x) {
    return sqrt (x);
  }
  static mag_type eps () {
    return epsilon ();
  }
  static mag_type sfmin () {
    return ArithTraits<mag_type>::sfmin ();
  }
  static int base () {
    return ArithTraits<mag_type>::base ();
  }
  static mag_type prec () {
    return ArithTraits<mag_type>::prec ();
  }
  static int t () {
    return ArithTraits<mag_type>::t ();
  }
  static mag_type rnd () {
    return ArithTraits<mag_type>::one ();
  }
  static int emin () {
    return ArithTraits<mag_type>::emin ();
  }
  static mag_type rmin () {
    return ArithTraits<mag_type>::rmin ();
  }
  static int emax () {
    return ArithTraits<mag_type>::emax ();
  }
  static mag_type rmax () {
    return ArithTraits<mag_type>::rmax ();
  }
};



template<>
class ArithTraits<double> {
public:
  typedef double mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const double x) {
#ifdef __CUDACC__
    return isinf (x);
#else
    return std::isinf (x);
#endif // __CUDACC__
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const double x) {
#ifdef __CUDACC__
    return isnan (x);
#else
    return std::isnan (x);
#endif // __CUDACC__
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const double x) {
    return ::fabs (x);
  }
  static KOKKOS_DEVICE_FUNCTION double zero () {
    return 0.0;
  }
  static KOKKOS_DEVICE_FUNCTION double one () {
    return 1.0;
  }
  static KOKKOS_DEVICE_FUNCTION double min () {
    return DBL_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION double max () {
    return DBL_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const double x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const double) {
    return 0.0;
  }
  static KOKKOS_DEVICE_FUNCTION double conj (const double x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION double pow (const double x, const double y) {
    return ::pow (x, y);
  }
  static KOKKOS_DEVICE_FUNCTION double sqrt (const double x) {
    return ::sqrt (x);
  }
  static KOKKOS_DEVICE_FUNCTION double log (const double x) {
    return ::log (x);
  }
  static KOKKOS_DEVICE_FUNCTION double log10 (const double x) {
    return ::log10 (x);
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return DBL_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef double magnitudeType;
  typedef float halfPrecision;
  typedef double doublePrecision; // CUDA doesn't support long double, unfortunately

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const double x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION double conjugate (const double x) {
    return conj (x);
  }
  //static KOKKOS_DEVICE_FUNCTION T nan ();
  //static KOKKOS_DEVICE_FUNCTION bool isnaninf (const double);
  static std::string name () {
    return "double";
  }
  static KOKKOS_DEVICE_FUNCTION double squareroot (const double x) {
    return sqrt (x);
  }
  static KOKKOS_DEVICE_FUNCTION double nan () {
#ifdef __CUDA_ARCH__
    return nan ();
#else
    // http://pubs.opengroup.org/onlinepubs/009696899/functions/nan.html
    return strtod ("NAN()", (char**) NULL);
#endif // __CUDA_ARCH__
  }
  static KOKKOS_DEVICE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_DEVICE_FUNCTION mag_type sfmin () {
    return DBL_MIN; // ???
  }
  static KOKKOS_DEVICE_FUNCTION int base () {
    return FLT_RADIX; // same for float as for double
  }
  static KOKKOS_DEVICE_FUNCTION mag_type prec () {
    return eps () * static_cast<mag_type> (base ());
  }
  static KOKKOS_DEVICE_FUNCTION int t () {
    return DBL_MANT_DIG;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type rnd () {
    return 1.0;
  }
  static KOKKOS_DEVICE_FUNCTION int emin () {
    return DBL_MIN_EXP;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type rmin () {
    return DBL_MIN; // ??? // should be base^(emin-1)
  }
  static KOKKOS_DEVICE_FUNCTION int emax () {
    return DBL_MAX_EXP;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type rmax () {
    return DBL_MAX; // ??? // should be (base^emax)*(1-eps)
  }
};



template<>
class ArithTraits<long double> {
public:
  typedef long double val_type;
  typedef long double mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static bool isInf (const val_type& x) {
#ifdef __CUDACC__
    return isinf (x);
#else
    return std::isinf (x);
#endif // __CUDACC__
  }
  static bool isNan (const val_type& x) {
#ifdef __CUDACC__
    return isnan (x);
#else
    return std::isnan (x);
#endif // __CUDACC__
  }
  static mag_type abs (const val_type& x) {
    return ::fabs (x);
  }
  static val_type zero () {
    return 0.0;
  }
  static val_type one () {
    return 1.0;
  }
  static val_type min () {
    return LDBL_MIN;
  }
  static val_type max () {
    return LDBL_MAX;
  }
  static mag_type real (const val_type& x) {
    return x;
  }
  static mag_type imag (const val_type&) {
    return zero ();
  }
  static val_type conj (const val_type& x) {
    return x;
  }
  static val_type pow (const val_type& x, const val_type& y) {
    return ::pow (x, y);
  }
  static val_type sqrt (const val_type& x) {
    return ::sqrt (x);
  }
  static val_type log (const val_type& x) {
    return ::log (x);
  }
  static val_type log10 (const val_type& x) {
    return ::log10 (x);
  }
  static mag_type epsilon () {
    return LDBL_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef val_type magnitudeType;
  typedef double halfPrecision;
  typedef val_type doublePrecision; // long double is the most you get, alas

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static mag_type magnitude (const val_type& x) {
    return abs (x);
  }
  static val_type conjugate (const val_type& x) {
    return conj (x);
  }
  //static KOKKOS_DEVICE_FUNCTION T nan ();
  //static KOKKOS_DEVICE_FUNCTION bool isnaninf (const double);
  static std::string name () {
    return "long double";
  }
  static val_type squareroot (const val_type& x) {
    return sqrt (x);
  }
  static val_type nan () {
    return strtold ("NAN()", (char**) NULL);
  }
  static mag_type eps () {
    return epsilon ();
  }
  static mag_type sfmin () {
    return LDBL_MIN; // ???
  }
  static int base () {
    return FLT_RADIX; // same for float as for double or long double
  }
  static mag_type prec () {
    return eps () * static_cast<mag_type> (base ());
  }
  static int t () {
    return LDBL_MANT_DIG;
  }
  static mag_type rnd () {
    return one ();
  }
  static int emin () {
    return LDBL_MIN_EXP;
  }
  static mag_type rmin () {
    return LDBL_MIN;
  }
  static int emax () {
    return LDBL_MAX_EXP;
  }
  static mag_type rmax () {
    return LDBL_MAX;
  }
};




// Interestingly enough, char and int8_t are different types, but
// signed char and int8_t are the same (on my system).  This means we
// need a separate specialization for char, but not for signed char or
// unsigned char.  Note that the C(++) standard does not specify
// whether char is signed or unsigned (!).
template<>
class ArithTraits<char> {
public:
  typedef char mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const char x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const char x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const char x) {
    // This may trigger a compiler warning if char is unsigned.
    // On most platforms I encounter, char is signed.
    return x >= 0 ? x : -x;
  }
  static KOKKOS_DEVICE_FUNCTION char zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION char one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION char min () {
    return CHAR_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION char max () {
    return CHAR_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const char x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const char x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION char conj (const char x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION char pow (const char x, const char y) {
    // if (y == 0) {
    //   if (x == 0) {
    //     // It's not entirely clear what to return if x and y are both
    //     // zero.  In the case of floating-point numbers, 0^0 is NaN.
    //     // Here, though, I think it's safe to return 0.
    //     return 0;
    //   }
    //   else {
    //     return 1;
    //   }
    // } else {
    //   char z = x; // skip the first iteration, since we know its result
    //   for (char k = 1; k < y; ++k) {
    //     z *= x;
    //   }
    //   return z;
    // }

    // We don't know if char is signed everywhere, but it generally is.
    return intPowSigned<char> (x, y);
  }
  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION char sqrt (const char x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<char> ( ::sqrt (static_cast<float> (abs (x))));
  }

  static KOKKOS_DEVICE_FUNCTION char log (const char x) {
    return static_cast<char> ( ::log (static_cast<float> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION char log10 (const char x) {
    return static_cast<char> ( ::log10 (static_cast<float> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef char magnitudeType;
  typedef char halfPrecision;
  typedef char doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const char x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION char conjugate (const char x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const char) {
    return false;
  }
  static std::string name () {
    return "char";
  }
  static KOKKOS_DEVICE_FUNCTION char squareroot (const char x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<int8_t> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef int8_t magnitudeType;
  typedef int8_t halfPrecision;
  typedef int8_t doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const int8_t x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION int8_t conjugate (const int8_t x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const int8_t) {
    return false;
  }
  static std::string name () {
    return "int8_t";
  }
  static KOKKOS_DEVICE_FUNCTION int8_t squareroot (const int8_t x) {
    return sqrt (x);
  }


  typedef int8_t mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const int8_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const int8_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const int8_t x) {
    // std::abs appears to work with CUDA 5.5 at least, but I'll use
    // the ternary expression for maximum generality.  Note that this
    // expression does not necessarily obey the rules for fabs() with
    // NaN input, so it should not be used for floating-point types.
    // It's perfectly fine for signed integer types, though.
    return x >= 0 ? x : -x;
  }
  static KOKKOS_DEVICE_FUNCTION int8_t zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int8_t one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION int8_t min () {
    // Macros like INT8_MIN work with CUDA, but
    // std::numeric_limits<int8_t>::min() does not, because it is not
    // marked as a __device__ function.
    return INT8_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION int8_t max () {
    return INT8_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const int8_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const int8_t x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int8_t conj (const int8_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION int8_t pow (const int8_t x, const int8_t y) {
    return intPowSigned<int8_t> (x, y);
  }
  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION int8_t sqrt (const int8_t x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<int8_t> ( ::sqrt (static_cast<float> (abs (x))));
  }

  static KOKKOS_DEVICE_FUNCTION int8_t log (const int8_t x) {
    return static_cast<int8_t> ( ::log (static_cast<float> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION int8_t log10 (const int8_t x) {
    return static_cast<int8_t> ( ::log10 (static_cast<float> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};



template<>
class ArithTraits<uint8_t> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef uint8_t magnitudeType;
  typedef uint8_t halfPrecision;
  typedef uint8_t doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const uint8_t x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t conjugate (const uint8_t x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const uint8_t) {
    return false;
  }
  static std::string name () {
    return "uint8_t";
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t squareroot (const uint8_t x) {
    return sqrt (x);
  }


  typedef uint8_t mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const uint8_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const uint8_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const uint8_t x) {
    return x; // it's unsigned, so it's positive
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t min () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t max () {
    return UINT8_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const uint8_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const uint8_t x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t conj (const uint8_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t pow (const uint8_t x, const uint8_t y) {
    return intPowUnsigned<uint8_t> (x, y);
  }

  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION uint8_t sqrt (const uint8_t x) {
    // It's not clear what to return if x is negative.  We could throw
    // an exception, but that won't work on the GPU.  If x were a
    // floating-point value, we could return NaN, but we can't do that
    // with integers.  So, we just take the absolute value and hope
    // for the best.

    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<uint8_t> ( ::sqrt (static_cast<float> (x)));
  }

  static KOKKOS_DEVICE_FUNCTION uint8_t log (const uint8_t x) {
    return static_cast<uint8_t> ( ::log (static_cast<float> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION uint8_t log10 (const uint8_t x) {
    return static_cast<uint8_t> ( ::log10 (static_cast<float> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};


template<>
class ArithTraits<int16_t> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef int16_t magnitudeType;
  typedef int16_t halfPrecision;
  typedef int16_t doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const int16_t x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION int16_t conjugate (const int16_t x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const int16_t) {
    return false;
  }
  static std::string name () {
    return "int16_t";
  }
  static KOKKOS_DEVICE_FUNCTION int16_t squareroot (const int16_t x) {
    return sqrt (x);
  }


  typedef int16_t mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const int16_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const int16_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const int16_t x) {
    // std::abs appears to work with CUDA 5.5 at least, but I'll use
    // the ternary expression for maximum generality.  Note that this
    // expression does not necessarily obey the rules for fabs() with
    // NaN input, so it should not be used for floating-point types.
    // It's perfectly fine for signed integer types, though.
    return x >= 0 ? x : -x;
  }
  static KOKKOS_DEVICE_FUNCTION int16_t zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int16_t one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION int16_t min () {
    // Macros like INT16_MIN work with CUDA, but
    // std::numeric_limits<int16_t>::min() does not, because it is not
    // marked as a __device__ function.
    return INT16_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION int16_t max () {
    return INT16_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const int16_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const int16_t x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int16_t conj (const int16_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION int16_t pow (const int16_t x, const int16_t y) {
    return intPowSigned<int16_t> (x, y);
  }
  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION int16_t sqrt (const int16_t x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<int16_t> ( ::sqrt (static_cast<float> (abs (x))));
  }

  static KOKKOS_DEVICE_FUNCTION int16_t log (const int16_t x) {
    return static_cast<int16_t> ( ::log (static_cast<float> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION int16_t log10 (const int16_t x) {
    return static_cast<int16_t> ( ::log10 (static_cast<float> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};



template<>
class ArithTraits<uint16_t> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef uint16_t magnitudeType;
  typedef uint16_t halfPrecision;
  typedef uint16_t doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const uint16_t x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t conjugate (const uint16_t x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const uint16_t) {
    return false;
  }
  static std::string name () {
    return "uint16_t";
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t squareroot (const uint16_t x) {
    return sqrt (x);
  }


  typedef uint16_t mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const uint16_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const uint16_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const uint16_t x) {
    return x; // it's unsigned, so it's positive
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t min () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t max () {
    return UINT16_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const uint16_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const uint16_t x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t conj (const uint16_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t pow (const uint16_t x, const uint16_t y) {
    return intPowUnsigned<uint16_t> (x, y);
  }

  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION uint16_t sqrt (const uint16_t x) {
    // It's not clear what to return if x is negative.  We could throw
    // an exception, but that won't work on the GPU.  If x were a
    // floating-point value, we could return NaN, but we can't do that
    // with integers.  So, we just take the absolute value and hope
    // for the best.

    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<uint16_t> ( ::sqrt (static_cast<float> (x)));
  }

  static KOKKOS_DEVICE_FUNCTION uint16_t log (const uint16_t x) {
    return static_cast<uint16_t> ( ::log (static_cast<float> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION uint16_t log10 (const uint16_t x) {
    return static_cast<uint16_t> ( ::log10 (static_cast<float> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};


template<>
class ArithTraits<int32_t> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef int32_t magnitudeType;
  typedef int32_t halfPrecision;
  typedef int32_t doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const int32_t x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION int32_t conjugate (const int32_t x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const int32_t) {
    return false;
  }
  static std::string name () {
    return "int32_t";
  }
  static KOKKOS_DEVICE_FUNCTION int32_t squareroot (const int32_t x) {
    return sqrt (x);
  }

  typedef int32_t mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const int32_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const int32_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const int32_t x) {
    // std::abs appears to work with CUDA 5.5 at least, but I'll use
    // the ternary expression for maximum generality.  Note that this
    // expression does not necessarily obey the rules for fabs() with
    // NaN input, so it should not be used for floating-point types.
    // It's perfectly fine for signed integer types, though.
    return x >= 0 ? x : -x;
  }
  static KOKKOS_DEVICE_FUNCTION int32_t zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int32_t one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION int32_t min () {
    // Macros like INT32_MIN work with CUDA, but
    // std::numeric_limits<int32_t>::min() does not, because it is not
    // marked as a __device__ function.
    return INT32_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION int32_t max () {
    return INT32_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const int32_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const int32_t x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int32_t conj (const int32_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION int32_t pow (const int32_t x, const int32_t y) {
    return intPowSigned<int32_t> (x, y);
  }
  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION int32_t sqrt (const int32_t x) {
    // if (x == 0) {
    //   return 0;
    // }
    // else {
    //   const int32_t x_abs = abs (x);
    //   // Square root of 2^{2k} is 2^k.
    //   const int32_t k = highestBit (x_abs); // floor(log2(abs(x)))
    //   return 1 << (k/2);
    // }

    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<int32_t> ( ::sqrt (static_cast<double> (abs (x))));
  }

  static KOKKOS_DEVICE_FUNCTION int32_t log (const int32_t x) {
    return static_cast<int32_t> ( ::log (static_cast<double> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION int32_t log10 (const int32_t x) {
    return static_cast<int32_t> ( ::log10 (static_cast<double> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};


template<>
class ArithTraits<uint32_t> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef uint32_t magnitudeType;
  typedef uint32_t halfPrecision;
  typedef uint32_t doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const uint32_t x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t conjugate (const uint32_t x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const uint32_t) {
    return false;
  }
  static std::string name () {
    return "uint32_t";
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t squareroot (const uint32_t x) {
    return sqrt (x);
  }

  typedef uint32_t mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const uint32_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const uint32_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const uint32_t x) {
    return x; // it's unsigned, so it's positive
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t min () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t max () {
    return UINT32_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const uint32_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const uint32_t x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t conj (const uint32_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t pow (const uint32_t x, const uint32_t y) {
    return intPowUnsigned<uint32_t> (x, y);
  }

  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION uint32_t sqrt (const uint32_t x) {
    // It's not clear what to return if x is negative.  We could throw
    // an exception, but that won't work on the GPU.  If x were a
    // floating-point value, we could return NaN, but we can't do that
    // with integers.  So, we just take the absolute value and hope
    // for the best.

    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<uint32_t> ( ::sqrt (static_cast<double> (x)));
  }

  static KOKKOS_DEVICE_FUNCTION uint32_t log (const uint32_t x) {
    return static_cast<uint32_t> ( ::log (static_cast<double> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION uint32_t log10 (const uint32_t x) {
    return static_cast<uint32_t> ( ::log10 (static_cast<double> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};


template<>
class ArithTraits<int64_t> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef int64_t magnitudeType;
  typedef int64_t halfPrecision;
  typedef int64_t doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const int64_t x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION int64_t conjugate (const int64_t x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const int64_t) {
    return false;
  }
  static std::string name () {
    return "int64_t";
  }
  static KOKKOS_DEVICE_FUNCTION int64_t squareroot (const int64_t x) {
    return sqrt (x);
  }

  typedef int64_t mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const int64_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const int64_t x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const int64_t x) {
    return x >= 0 ? x : -x;
  }
  static KOKKOS_DEVICE_FUNCTION int64_t zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int64_t one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION int64_t min () {
    return INT64_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION int64_t max () {
    return INT64_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const int64_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const int64_t x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION int64_t conj (const int64_t x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION int64_t pow (const int64_t x, const int64_t y) {
    return intPowSigned<int64_t> (x, y);
  }
  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION int64_t sqrt (const int64_t x) {
    // Casting from int64_t to double does result in a loss of
    // accuracy.  However, it gives us a good first approximation.
    // For very large numbers, we may lose some significand bits, but
    // will always get within a factor of two (assuming correct
    // rounding) of the exact double-precision number.  We could then
    // binary search between half the result and twice the result
    // (assuming the latter is <= INT64_MAX, which it has to be, so we
    // don't have to check) to ensure correctness.  It actually should
    // suffice to check numbers within 1 of the result.

    const int64_t approx = static_cast<int64_t> ( ::sqrt (static_cast<double> (abs (x))));
    // const int64_t approx_squared = approx * approx;
    // const int64_t approx_plus_one_squared = (approx+1) * (approx+1);
    // const int64_t approx_minus_one_squared = (approx-1) * (approx-1);

    // if (approx_squared > x) {
    //   return approx - 1;
    // } else if (approx_squared < x) {
    //   if (approx_plus_one_squared < x) {
    //     return approx + 1;
    //   } else {
    //     return approx;
    //   }
    // } else {
    //   return approx; // exactly right
    // }

    return approx;
  }

  static KOKKOS_DEVICE_FUNCTION int64_t log (const int64_t x) {
    return static_cast<int64_t> ( ::log (static_cast<double> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION int64_t log10 (const int64_t x) {
    return static_cast<int64_t> ( ::log10 (static_cast<double> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};


// With CUDA 5.5, long and int64_t are different, even though
// sizeof(long) == 8.  This could be because long and long long are
// two different types (not aliases of one another), even though they
// might have the same size.
template<>
class ArithTraits<long> {
public:
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef long magnitudeType;
  typedef long halfPrecision;
  typedef long doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const long x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION long conjugate (const long x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const long) {
    return false;
  }
  static std::string name () {
    return "long";
  }
  static KOKKOS_DEVICE_FUNCTION long squareroot (const long x) {
    return sqrt (x);
  }

  typedef long mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const long x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const long x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const long x) {
    return x >= 0 ? x : -x;
  }
  static KOKKOS_DEVICE_FUNCTION long zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION long one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION long min () {
    return LONG_MIN;
  }
  static KOKKOS_DEVICE_FUNCTION long max () {
    return LONG_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const long x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const long x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION long conj (const long x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION long pow (const long x, const long y) {
    return intPowSigned<long> (x, y);
  }
  static KOKKOS_DEVICE_FUNCTION long sqrt (const long x) {
    return static_cast<long> ( ::sqrt (static_cast<double> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION long log (const long x) {
    return static_cast<long> ( ::log (static_cast<double> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION long log10 (const long x) {
    return static_cast<long> ( ::log10 (static_cast<double> (abs (x))));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }
};


template<>
class ArithTraits<unsigned long> {
public:
  typedef unsigned long val_type;
  typedef unsigned long mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const val_type x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const val_type x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const val_type x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION val_type min () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION val_type max () {
    return ULONG_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const val_type x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION val_type pow (const val_type x, const val_type y) {
    return intPowUnsigned<val_type> (x, y);
  }
  static KOKKOS_DEVICE_FUNCTION val_type sqrt (const val_type x) {
    return static_cast<val_type> ( ::sqrt (static_cast<double> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION val_type log (const val_type x) {
    return static_cast<long> ( ::log (static_cast<double> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<long> ( ::log10 (static_cast<double> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef val_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "unsigned long";
  }
  static KOKKOS_DEVICE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};



template<>
class ArithTraits<uint64_t> {
public:
  typedef uint64_t val_type;
  typedef uint64_t mag_type;
  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static KOKKOS_DEVICE_FUNCTION bool isInf (const val_type x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION bool isNan (const val_type x) {
    return false;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type abs (const val_type x) {
    return x; // unsigned integers are always nonnegative
  }
  static KOKKOS_DEVICE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_DEVICE_FUNCTION val_type min () {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION val_type max () {
    return UINT64_MAX;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION mag_type imag (const val_type x) {
    return 0;
  }
  static KOKKOS_DEVICE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_DEVICE_FUNCTION val_type pow (const val_type x, const val_type y) {
    return intPowUnsigned<val_type> (x, y);
  }
  //! Integer square root returns a lower bound.
  static KOKKOS_DEVICE_FUNCTION val_type sqrt (const val_type x) {
    // Casting from uint64_t to to double does result in a loss of
    // accuracy.  However, it gives us a good first approximation.
    // For very large numbers, we may lose some significand bits, but
    // will always get within a factor of two (assuming correct
    // rounding) of the exact double-precision number.  We could then
    // binary search between half the result and twice the result
    // (assuming the latter is <= UINT64_MAX, which it has to be, so
    // we don't have to check) to ensure correctness.  It actually
    // should suffice to check numbers within 1 of the result.
    return static_cast<val_type> ( ::sqrt (static_cast<double> (x)));
  }

  static KOKKOS_DEVICE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<double> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<double> (x)));
  }
  static KOKKOS_DEVICE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef val_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_DEVICE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_DEVICE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_DEVICE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "uint64_t";
  }
  static KOKKOS_DEVICE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};

} // namespace Details
} // namespace Kokkos

#endif // KOKKOS_ARITHTRAITS_HPP
