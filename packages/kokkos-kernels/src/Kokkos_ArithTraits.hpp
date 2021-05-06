/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_ARITHTRAITS_HPP
#define KOKKOS_ARITHTRAITS_HPP

/// \file Kokkos_ArithTraits.hpp
/// \brief Declaration and definition of Kokkos::Details::ArithTraits

#include <KokkosKernels_config.h>
#include <Kokkos_Complex.hpp>
#include <KokkosKernels_Half.hpp>
#include <Kokkos_Macros.hpp>

#ifdef HAVE_KOKKOSKERNELS_QUADMATH
#  include <quadmath.h>
#endif // HAVE_KOKKOSKERNELS_QUADMATH

#include <cfloat>
#include <climits>
#include <cmath>
#include <complex> // std::complex
#include <limits> // std::numeric_limits
#ifdef __CUDACC__
#  include <math_constants.h>
#endif

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
KOKKOS_FORCEINLINE_FUNCTION IntType
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


// Warning free abs function for types where we don't know whether they are signed (like char)
template<class T, bool is_signed = std::numeric_limits<T>::is_signed >
struct integer_abs {
  static
  KOKKOS_INLINE_FUNCTION
  T abs(const T& val);
};

template<class T>
struct integer_abs<T,true> {
  static
  KOKKOS_INLINE_FUNCTION
  T abs(const T& x) {
    return x<0? -x:x;
  }
};

template<class T>
struct integer_abs<T,false> {
  static
  KOKKOS_INLINE_FUNCTION
  T abs(const T& x) {
    return x;
  }
};



/// \fn intPowSigned
/// \tparam IntType A built-in signed integer type.
/// \brief Compute x raised to the power y.
///
/// If the arguments are invalid (e.g., if x and y are both zero), the
/// result of this function is undefined.  However, this function will
/// not throw an exception in that case.
template<class IntType>
KOKKOS_FORCEINLINE_FUNCTION 
typename std::enable_if<std::numeric_limits<IntType>::is_signed,IntType>::type
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
  }
  return intPowImpl<IntType> (x, y);
}
template<class IntType>
KOKKOS_FORCEINLINE_FUNCTION 
typename std::enable_if<!std::numeric_limits<IntType>::is_signed,IntType>::type
intPowSigned (const IntType x, const IntType y)
{
  // It's not entirely clear what to return if x and y are both zero.
  // In the case of floating-point numbers, 0^0 is NaN.  Here, though,
  // I think it's safe to return 0.
  if (x == 0) {
    return 0;
  } else if (y == 0) {
    return 1;
  }
  return intPowImpl<IntType> (x, y);
}

/// \fn intPowUnsigned
/// \tparam IntType A built-in unsigned integer type.
/// \brief Compute x raised to the power y.
///
/// If the arguments are invalid (e.g., if x and y are both zero), the
/// result of this function is undefined.  However, this function will
/// not throw an exception in that case.
template<class IntType>
KOKKOS_FORCEINLINE_FUNCTION IntType
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

} // namespace (anonymous)

namespace Kokkos {
namespace Details {

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
template<class T>
class ArithTraits {
public:
  /// \brief A type that acts like T and works with Kokkos.
  ///
  /// This is usually just an alias for T.  However, some types T do
  /// not work well with Kokkos.  In that case, we use a mostly
  /// equivalent type here.  For example, ArithTraits<std::complex<R>
  /// >::val_type is Kokkos::complex<R>.
  typedef T val_type;
  /// \brief The type of the magnitude (absolute value) of T.
  ///
  /// We define this as the type returned by abs() in this class.  If
  /// T is real (not complex), then \c val_type and \c mag_type are
  /// usually the same.  If T is <tt>std::complex<R></tt> for some R,
  /// then R and \c mag_type are usually the same.
  typedef T mag_type;

  //! Whether ArithTraits has a specialization for T.
  static const bool is_specialized = false;
  //! Whether T is a signed type (has negative values).
  static const bool is_signed = false;
  //! Whether T is an integer type.
  static const bool is_integer = false;
  /// \brief Whether T "uses exact representations."
  ///
  /// The opposite of is_exact is "is approximate," that is, "may
  /// commit rounding error."
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
  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const T& x);

  /// \brief Whether x is NaN (not a number).
  ///
  /// This can only be true for floating-point types T that support
  /// NaN.  If T is a complex type, we say that a T instance x is NaN
  /// if and only if <tt>isNan(real(x)) || isNan(imag(x))</tt>.
  ///
  /// Unfortunately we can't call this "isnan" (the equivalent C99
  /// function), because CUDA appears to implement that function using
  /// a macro, rather than using a function (as C++11 requires).
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const T& x);

  //! The absolute value (magnitude) of x.
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const T& x);

  //! The zero value of T; the arithmetic identity.
  static KOKKOS_FORCEINLINE_FUNCTION T zero ();

  //! The one value of T; the multiplicative identity.
  static KOKKOS_FORCEINLINE_FUNCTION T one ();

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
  static KOKKOS_FORCEINLINE_FUNCTION T infinity();

  /// \brief The minimum possible value of T.
  ///
  /// If T is a real floating-point type, then this is the minimum
  /// <i>positive</i> value, as with std::numeric_limits<T>::min().
  static KOKKOS_FORCEINLINE_FUNCTION T min ();

  //! The maximum possible value of T.
  static KOKKOS_FORCEINLINE_FUNCTION T max ();

  /// \brief The real part of x.
  ///
  /// If \c is_complex is false, then this just returns x.
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const T& x);

  /// \brief The imaginary part of x.
  ///
  /// If \c is_complex is false, then this just returns zero().
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const T&);

  /// \brief The complex conjugate of x.
  ///
  /// If \c is_complex is false, then this just returns x.
  static KOKKOS_FORCEINLINE_FUNCTION T conj (const T&);

  //! x raised to the power y.
  static KOKKOS_FORCEINLINE_FUNCTION T pow (const T& x, const T& y);

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
  static KOKKOS_FORCEINLINE_FUNCTION T sqrt (const T& x);

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
  static KOKKOS_FORCEINLINE_FUNCTION T cbrt (const T& x);  

  /// \brief The natural (base e) exponential function of x.
  ///
  /// If T is an integer type, this is the floor of the exponential 
  /// function.  If T is a complex-valued type, then this method 
  /// returns \f$e^{x+iy} = e^x ( cos(y) + i sin(y) )\f$.
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T exp (const T& x);

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
  static KOKKOS_FORCEINLINE_FUNCTION T log (const T& x);

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
  static KOKKOS_FORCEINLINE_FUNCTION T log10 (const T& x);

  /// Trigonometric and hyperbolic functions are not available 
  /// for integer types. This is because asin(sin(x)) is not x
  /// when x is integer with a rounding error.
  ///
  ///  KJ: log, exp also has this problem. We probably need to 
  ///      disable them for integer types instead of providing 
  ///      functionality with floor.

  /// \brief The sin function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T sin (const T& x);

  /// \brief The cos function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T cos (const T& x);

  /// \brief The tan function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T tan (const T& x);
  
  /// \brief The sin hyperbolic function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T sinh (const T& x);

  /// \brief The cos hyperbolic function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T cosh (const T& x);

  /// \brief The tan hyperbolic function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T tanh (const T& x);
  
  /// \brief The asin function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T asin (const T& x);

  /// \brief The acos function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T acos (const T& x);

  /// \brief The atan function of x
  ///
  static KOKKOS_FORCEINLINE_FUNCTION T atan (const T& x);

  /// \brief Return a silent NaN, if appropriate for T.
  ///
  /// If T does <i>not</i> implement a silent NaN, the return value is
  /// undefined, but calling this method is still allowed.
  static KOKKOS_FORCEINLINE_FUNCTION T nan ();

  /// \brief Machine epsilon.
  ///
  /// If T is an integer type (std::numeric_traits<T>::is_exact is
  /// true), then epsilon() returns 0.  Otherwise, if T is a
  /// floating-point type, it returns machine epsilon that T.
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon ();

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
  static KOKKOS_FORCEINLINE_FUNCTION mag_type eps ();

  //! Return safe minimum (sfmin), such that 1/sfmin does not overflow.
  static KOKKOS_FORCEINLINE_FUNCTION mag_type sfmin ();

  //! Return the base of the scalar type T.
  static KOKKOS_FORCEINLINE_FUNCTION int base ();

  //! Return <tt>eps*base</tt>.
  static KOKKOS_FORCEINLINE_FUNCTION mag_type prec ();

  //! Returns the number of (base) digits in the significand.
  static KOKKOS_FORCEINLINE_FUNCTION int t ();

  //! 1.0 when rounding occurs in addition, else 0.0.
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rnd ();

  //! Returns the minimum exponent before (gradual) underflow.
  static KOKKOS_FORCEINLINE_FUNCTION int emin ();

  //! Returns the underflow threshold: <tt>base^(emin-1)</tt>
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmin ();

  //! Returns the largest exponent before overflow.
  static KOKKOS_FORCEINLINE_FUNCTION int emax ();

  //! Overflow theshold: <tt>(base^emax)*(1-eps)</tt>
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmax ();

  //! Same as abs(); return the magnitude of x.
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const T& x);

  //! Same as conj(); return the complex conjugate of x.
  static KOKKOS_FORCEINLINE_FUNCTION T conjugate (const T& x);

  /// \brief Whether x is (silent) NaN or Inf.
  ///
  /// This is the same as <tt>isNan(x) || isInf(x)</tt>.
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const T& x);

  /// \brief The string name of T.
  ///
  /// Note that this is not a device function.
  static std::string name ();

  //! Same as sqrt(x); the square root of x.
  static KOKKOS_FORCEINLINE_FUNCTION T squareroot (const T& x);
  //@}
};

// Since Kokkos::Experimental::half_t falls back to float, only define
// ArithTraits if half_t is a backend specialization
#if defined(KOKKOS_HALF_T_IS_FLOAT) &&\
    !KOKKOS_HALF_T_IS_FLOAT
template <>
class ArithTraits<Kokkos::Experimental::half_t> {
public:
  typedef Kokkos::Experimental::half_t val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static constexpr bool has_infinity = true;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return Kokkos::Experimental::cast_to_half(HUGE_VALF); }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type x) {
    #ifndef __CUDA_ARCH__
    using std::isinf;
    #endif
    return isinf (Kokkos::Experimental::cast_from_half<float>(x));
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type x) {
    #ifndef __CUDA_ARCH__
    using std::isnan;
    #endif
    return isnan(Kokkos::Experimental::cast_from_half<float>(x));
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return Kokkos::Experimental::cast_to_half(fabs(Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return Kokkos::Experimental::cast_to_half(0.0F);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return Kokkos::Experimental::cast_to_half(1.0F);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return Kokkos::Experimental::cast_to_half(-KOKKOSKERNELS_IMPL_FP16_MAX);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return Kokkos::Experimental::cast_to_half(KOKKOSKERNELS_IMPL_FP16_MAX);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type) {
    return Kokkos::Experimental::cast_to_half(0.0F);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const val_type y) {
    return Kokkos::Experimental::cast_to_half(::pow(Kokkos::Experimental::cast_from_half<float>(x),
                 Kokkos::Experimental::cast_from_half<float>(y)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::sqrt (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::cbrt (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::exp (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::log (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::log10 (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::sin (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::cos (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::tan (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::sinh (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::cosh (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::tanh (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::asin (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::acos (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
    return Kokkos::Experimental::cast_to_half(::atan (Kokkos::Experimental::cast_from_half<float>(x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    //return ::pow(2, -KOKKOSKERNELS_IMPL_FP16_SIGNIFICAND_BITS);
    return Kokkos::Experimental::cast_to_half(KOKKOSKERNELS_IMPL_FP16_EPSILON);
  }
  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  // C++ doesn't have a standard "half-float" type.
  typedef val_type halfPrecision;
  typedef double doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type x) {
    return isNan (x) || isInf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs(x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj(x);
  }
  static std::string name () {
    return "half";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt(x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
#ifdef __CUDA_ARCH__
    return Kokkos::Experimental::cast_to_half(CUDART_NAN_F);
#else
    return Kokkos::Experimental::cast_to_half(std::numeric_limits<float>::quiet_NaN());
#endif // __CUDA_ARCH__
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type sfmin () {
    return Kokkos::Experimental::cast_to_half(KOKKOSKERNELS_IMPL_FP16_MIN);
  }
  static KOKKOS_FORCEINLINE_FUNCTION int base () {
    return KOKKOSKERNELS_IMPL_FP16_RADIX;
  }
  // Use float to allow running on both host and device
  static KOKKOS_FORCEINLINE_FUNCTION float prec () {
    float e = KOKKOSKERNELS_IMPL_FP16_EPSILON;
    float b = (float) base();
    float r = e * b;
    return r;
  }
  static KOKKOS_FORCEINLINE_FUNCTION int t () {
    return KOKKOSKERNELS_IMPL_FP16_MANT_DIG;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rnd () {
    return Kokkos::Experimental::cast_to_half(1.0);
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emin () {
    return KOKKOSKERNELS_IMPL_FP16_MIN_EXP;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmin () {
    return Kokkos::Experimental::cast_to_half(KOKKOSKERNELS_IMPL_FP16_MIN);
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emax () {
    return KOKKOSKERNELS_IMPL_FP16_MAX_EXP;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmax () {
    return Kokkos::Experimental::cast_to_half(KOKKOSKERNELS_IMPL_FP16_MAX);
  }
};
#endif // KOKKOS_HALF_T_IS_FLOAT && KOKKOS_ENABLE_CUDA_HALF

template<>
class ArithTraits<float> {
public:
  typedef float val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static constexpr bool has_infinity = true;
  static KOKKOS_FORCEINLINE_FUNCTION float infinity() { return HUGE_VALF; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const float x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isinf;
#elif KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    using sycl::isinf
#endif
    return isinf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const float x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isnan;
#elif KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    using sycl::isnan
#endif
    return isnan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const float x) {
    return ::fabs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float zero () {
    return 0.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION float one () {
    return 1.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION float min () {
    return -FLT_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION float max () {
    return FLT_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const float x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const float) {
    return 0.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION float conj (const float x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION float pow (const float x, const float y) {
    return ::pow (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float sqrt (const float x) {
    return ::sqrt (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float cbrt (const float x) {
    return ::cbrt (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float exp (const float x) {
    return ::exp (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float log (const float x) {
    return ::log (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float log10 (const float x) {
    return ::log10 (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float sin (const float x) {
    return ::sin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float cos (const float x) {
    return ::cos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float tan (const float x) {
    return ::tan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float sinh (const float x) {
    return ::sinh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float cosh (const float x) {
    return ::cosh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float tanh (const float x) {
    return ::tanh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float asin (const float x) {
    return ::asin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float acos (const float x) {
    return ::acos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float atan (const float x) {
    return ::atan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return FLT_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  // C++ doesn't have a standard "half-float" type.
  typedef float halfPrecision;
  typedef double doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const float x) {
    return isNan (x) || isInf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const float x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float conjugate (const float x) {
    return conj (x);
  }
  static std::string name () {
    return "float";
  }
  static KOKKOS_FORCEINLINE_FUNCTION float squareroot (const float x) {
    return sqrt (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION float nan () {
#if defined(__CUDA_ARCH__)
    return CUDART_NAN_F;
    //return nan (); //this returns 0???
#elif defined(__HIP_DEVICE_COMPILE__)
    return ::nanf("");
#else
    return std::numeric_limits<float>::quiet_NaN();
#endif // __CUDA_ARCH__
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type sfmin () {
    return FLT_MIN; // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION int base () {
    return FLT_RADIX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type prec () {
    return eps () * static_cast<mag_type> (base ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION int t () {
    return FLT_MANT_DIG;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rnd () {
    return 1.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emin () {
    return FLT_MIN_EXP;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmin () {
    return FLT_MIN; // ??? // should be base^(emin-1)
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emax () {
    return FLT_MAX_EXP;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmax () {
    return FLT_MAX; // ??? // should be (base^emax)*(1-eps)
  }
};

/// \brief Partial specialization for std::complex<RealFloatType>.
///
/// The C++ Standard Library (with C++03 at least) only allows
/// std::complex<RealFloatType> for RealFloatType = float, double, or
/// long double.
template<class RealFloatType>
class ArithTraits<std::complex<RealFloatType> > {
public:
  //! Kokkos internally replaces std::complex with Kokkos::complex.
  typedef ::Kokkos::complex<RealFloatType> val_type;
  typedef RealFloatType mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = true;

  static constexpr bool has_infinity = true;
  static std::complex<RealFloatType> infinity() {
    return std::complex<RealFloatType> (ArithTraits<mag_type>::infinity (), ArithTraits<mag_type>::infinity ());
  }

#ifdef KOKKOS_ENABLE_SYCL
  template <typename Dummy = RealFloatType>
  static bool isInf(const std::complex<Dummy>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isinf;
#elif KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    using sycl::isinf
#endif
    return isinf (real (x)) || isinf (imag (x));
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
    return isinf (real (x)) || isinf (imag (x));
  }
#endif
#ifdef KOKKOS_ENABLE_SYCL
  template <typename Dummy = RealFloatType>
  static bool isNan(const std::complex<Dummy>& x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isnan;
#elif KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    using sycl::isnan
#endif
    return isnan (real (x)) || isnan (imag (x));
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
    return isnan (real (x)) || isnan (imag (x));
  }
#endif
  static mag_type abs (const std::complex<RealFloatType>& x) {
    return std::abs (x);
  }
  static std::complex<RealFloatType> zero () {
    return std::complex<RealFloatType> (ArithTraits<mag_type>::zero (), ArithTraits<mag_type>::zero ());
  }
  static std::complex<RealFloatType> one () {
    return std::complex<RealFloatType> (ArithTraits<mag_type>::one (), ArithTraits<mag_type>::zero ());
  }
  static std::complex<RealFloatType> min () {
    return std::complex<RealFloatType> (ArithTraits<mag_type>::min (), ArithTraits<mag_type>::zero ());
  }
  static std::complex<RealFloatType> max () {
    return std::complex<RealFloatType> (ArithTraits<mag_type>::max (), ArithTraits<mag_type>::zero ());
  }
  static mag_type real (const std::complex<RealFloatType>& x) {
    return std::real (x);
  }
  static mag_type imag (const std::complex<RealFloatType>& x) {
    return std::imag (x);
  }
  static std::complex<RealFloatType> conj (const std::complex<RealFloatType>& x) {
    return std::conj (x);
  }
  static std::complex<RealFloatType>
  pow (const std::complex<RealFloatType>& x, const std::complex<RealFloatType>& y) {
    // Fix for some weird gcc 4.2.1 inaccuracy.
    if (y == one ()) {
      return x;
    } else if (y == one () + one ()) {
      return x * x;
    } else {
      return std::pow (x, y);
    }
  }
  static std::complex<RealFloatType>
  pow (const std::complex<RealFloatType>& x, const RealFloatType & y) {
    // Fix for some weird gcc 4.2.1 inaccuracy.
    if (y == ArithTraits<RealFloatType>::one ()) {
      return x;
    } else if (y == ArithTraits<RealFloatType>::one () + ArithTraits<RealFloatType>::one ()) {
      return x * x;
    } else {
      return std::pow (x, y);
    }
  }
  static std::complex<RealFloatType> sqrt (const std::complex<RealFloatType>& x) {
    return std::sqrt (x);
  }
  static std::complex<RealFloatType> cbrt (const std::complex<RealFloatType>& x) {
    return std::cbrt (x);
  }
  static std::complex<RealFloatType> exp (const std::complex<RealFloatType>& x) {
    return std::exp (x);
  }
  static std::complex<RealFloatType> log (const std::complex<RealFloatType>& x) {
    return std::log (x);
  }
  static std::complex<RealFloatType> log10 (const std::complex<RealFloatType>& x) {
    return std::log10 (x);
  }
  static std::complex<RealFloatType> sin (const std::complex<RealFloatType>& x) {
    return std::sin (x);
  }
  static std::complex<RealFloatType> cos (const std::complex<RealFloatType>& x) {
    return std::cos (x);
  }
  static std::complex<RealFloatType> tan (const std::complex<RealFloatType>& x) {
    return std::tan (x);
  }
  static std::complex<RealFloatType> sinh (const std::complex<RealFloatType>& x) {
    return std::sinh (x);
  }
  static std::complex<RealFloatType> cosh (const std::complex<RealFloatType>& x) {
    return std::cosh (x);
  }
  static std::complex<RealFloatType> tanh (const std::complex<RealFloatType>& x) {
    return std::tanh (x);
  }
  static std::complex<RealFloatType> asin (const std::complex<RealFloatType>& x) {
    return std::asin (x);
  }
  static std::complex<RealFloatType> acos (const std::complex<RealFloatType>& x) {
    return std::acos (x);
  }
  static std::complex<RealFloatType> atan (const std::complex<RealFloatType>& x) {
    return std::atan (x);
  }
  static std::complex<RealFloatType> nan () {
    const mag_type mag_nan = ArithTraits<mag_type>::nan ();
    return std::complex<RealFloatType> (mag_nan, mag_nan);
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
  static bool isnaninf (const std::complex<RealFloatType>& x) {
    return isNan (x) || isInf (x);
  }
  static mag_type magnitude (const std::complex<RealFloatType>& x) {
    return abs (x);
  }
  static std::complex<RealFloatType> conjugate (const std::complex<RealFloatType>& x) {
    return conj (x);
  }
  static std::string name () {
    return std::string ("std::complex<") + ArithTraits<mag_type>::name () + ">";
  }
  static std::complex<RealFloatType> squareroot (const std::complex<RealFloatType>& x) {
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
  typedef double val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static constexpr bool has_infinity = true;
  static KOKKOS_FORCEINLINE_FUNCTION double infinity() { return HUGE_VAL; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type x) {
    #ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isinf;
    #elif KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    using sycl::isinf;
    #endif
    return isinf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type x) {
    #ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::isnan;
    #elif KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL
    using sycl::isnan;
    #endif
    return isnan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return ::fabs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return -DBL_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return DBL_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type) {
    return 0.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const val_type y) {
    return ::pow (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    return ::sqrt (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return ::cbrt (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return ::exp (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return ::log (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return ::log10 (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
    return ::sin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
    return ::cos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
    return ::tan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
    return ::sinh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
    return ::cosh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
    return ::tanh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
    return ::asin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
    return ::acos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
    return ::atan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
#if defined(__CUDA_ARCH__)
    return CUDART_NAN;
    //return nan (); // this returns 0 ???
#elif defined(__HIP_DEVICE_COMPILE__)
    return ::nan("");
#else
    return std::numeric_limits<val_type>::quiet_NaN();
#endif // __CUDA_ARCH__
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return DBL_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef float halfPrecision;
#if defined(__CUDA_ARCH__)
  typedef double doublePrecision; // CUDA doesn't support long double, unfortunately
#elif defined(__HIP_DEVICE_COMPILE__)
  typedef double doublePrecision; // HIP does not support long double unfortunately
#else
  typedef long double doublePrecision;
#endif // __CUDA_ARCH__
  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static bool isnaninf (const val_type& x) {
    return isNan (x) || isInf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static std::string name () {
    return "double";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type sfmin () {
    return DBL_MIN; // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION int base () {
    return FLT_RADIX; // same for float as for double
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type prec () {
    return eps () * static_cast<mag_type> (base ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION int t () {
    return DBL_MANT_DIG;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rnd () {
    return 1.0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emin () {
    return DBL_MIN_EXP;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmin () {
    return DBL_MIN; // ??? // should be base^(emin-1)
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emax () {
    return DBL_MAX_EXP;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmax () {
    return DBL_MAX; // ??? // should be (base^emax)*(1-eps)
  }
};


// CUDA and HIP do not support long double in device functions,
// so none of the class methods in this specialization are marked
// as device functions.
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
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

  static constexpr bool has_infinity = true;
  static long double infinity() { return HUGE_VALL; }

  static bool isInf (const val_type& x) {
    using std::isinf;
    return isinf (x);
  }
  static bool isNan (const val_type& x) {
    using std::isnan;
    return isnan (x);
  }
  static mag_type abs (const val_type& x) {
    return ::fabsl (x);
  }
  static val_type zero () {
    return 0.0;
  }
  static val_type one () {
    return 1.0;
  }
  static val_type min () {
    return -LDBL_MAX;
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
  static val_type cbrt (const val_type& x) {
    return ::cbrtl (x);
  }
  static val_type exp (const val_type& x) {
    return ::exp (x);
  }
  static val_type log (const val_type& x) {
    return ::log (x);
  }
  static val_type log10 (const val_type& x) {
    return ::log10 (x);
  }
  static val_type sin (const val_type& x) {
    return ::sin (x);
  }
  static val_type cos (const val_type& x) {
    return ::cos (x);
  }
  static val_type tan (const val_type& x) {
    return ::tan (x);
  }
  static val_type sinh (const val_type& x) {
    return ::sinh (x);
  }
  static val_type cosh (const val_type& x) {
    return ::cosh (x);
  }
  static val_type tanh (const val_type& x) {
    return ::tanh (x);
  }
  static val_type asin (const val_type& x) {
    return ::asin (x);
  }
  static val_type acos (const val_type& x) {
    return ::acos (x);
  }
  static val_type atan (const val_type& x) {
    return ::atan (x);
  }
  static val_type nan () {
    return std::numeric_limits<val_type>::quiet_NaN();
  }
  static mag_type epsilon () {
    return LDBL_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef double halfPrecision;
  // It might be appropriate to use QD's qd_real here.
  // For now, long double is the most you get.
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
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
  static std::string name () {
    return "long double";
  }
  static val_type squareroot (const val_type& x) {
    return sqrt (x);
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
}; // long double specialization
#endif // KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST

#ifdef HAVE_KOKKOSKERNELS_QUADMATH

// CUDA does not support __float128 in device functions, so none of
// the class methods in this specialization are marked as device
// functions.
template<>
class ArithTraits<__float128> {
public:
  typedef __float128 val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static constexpr bool has_infinity = true;
  static __float128 infinity() { return 1.0q / 0.0q; }

  static bool isInf (const __float128 x) {
    return isinfq (x);
  }
  static bool isNan (const __float128 x) {
    return isnanq (x);
  }
  static mag_type abs (const __float128 x) {
    return fabsq (x);
  }
  static __float128 zero () {
    return 0.0;
  }
  static __float128 one () {
    return 1.0;
  }
  static __float128 min () {
    return FLT128_MIN;
  }
  static __float128 max () {
    return FLT128_MAX;
  }
  static mag_type real (const __float128 x) {
    return x;
  }
  static mag_type imag (const __float128 /* x */) {
    return 0.0;
  }
  static __float128 conj (const __float128 x) {
    return x;
  }
  static __float128 pow (const __float128 x, const __float128 y) {
    return powq (x, y);
  }
  static __float128 sqrt (const __float128 x) {
    return sqrtq (x);
  }
  static __float128 cbrt (const __float128 x) {
    return cbrtq (x);
  }
  static __float128 exp (const __float128 x) {
    return exp (x);
  }
  static __float128 log (const __float128 x) {
    return logq (x);
  }
  static __float128 log10 (const __float128 x) {
    return log10q (x);
  }
  static __float128 sin (const __float128 x) {
    return sinq (x);
  }
  static __float128 cos (const __float128 x) {
    return cosq (x);
  }
  static __float128 tan (const __float128 x) {
    return tanq (x);
  }
  static __float128 sinh (const __float128 x) {
    return sinhq (x);
  }
  static __float128 cosh (const __float128 x) {
    return coshq (x);
  }
  static __float128 tanh (const __float128 x) {
    return tanhq (x);
  }
  static __float128 asin (const __float128 x) {
    return asinq (x);
  }
  static __float128 acos (const __float128 x) {
    return acosq (x);
  }
  static __float128 atan (const __float128 x) {
    return atanq (x);
  }
  static mag_type epsilon () {
    return FLT128_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef double halfPrecision;
  // Unfortunately, we can't rely on a standard __float256 type.
  typedef __float128 doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static bool isnaninf (const __float128 x) {
    return isNan (x) || isInf (x);
  }
  static magnitudeType magnitude (const __float128 x) {
    return abs (x);
  }
  static __float128 conjugate (const __float128 x) {
    return conj (x);
  }
  static std::string name () {
    return "__float128";
  }
  static __float128 squareroot (const __float128 x) {
    return sqrt (x);
  }
  static __float128 nan () {
    return strtoflt128 ("NAN()", NULL); // ???
  }
  static mag_type eps () {
    return epsilon ();
  }
  static mag_type sfmin () {
    return FLT128_MIN; // ???
  }
  static int base () {
    return 2;
  }
  static mag_type prec () {
    return eps () * static_cast<mag_type> (base ());
  }
  static int t () {
    return FLT_MANT_DIG;
  }
  static mag_type rnd () {
    return 1.0;
  }
  static int emin () {
    return FLT128_MIN_EXP;
  }
  static mag_type rmin () {
    return FLT128_MIN; // ??? // should be base^(emin-1)
  }
  static int emax () {
    return FLT128_MAX_EXP;
  }
  static mag_type rmax () {
    return FLT128_MAX; // ??? // should be (base^emax)*(1-eps)
  }
};
#endif // HAVE_KOKKOSKERNELS_QUADMATH

template<>
class ArithTraits< ::Kokkos::complex<float> > {
public:
  typedef ::Kokkos::complex<float> val_type;
  typedef float mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = true;

  static constexpr bool has_infinity = true;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() {
    return val_type (ArithTraits<mag_type>::infinity (), ArithTraits<mag_type>::infinity ());
  }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type x) {
    return ArithTraits<mag_type>::isInf (x.real ()) ||
      ArithTraits<mag_type>::isInf (x.imag ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type x) {
    return ArithTraits<mag_type>::isNan (x.real ()) ||
      ArithTraits<mag_type>::isNan (x.imag ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return std::sqrt (::Kokkos::real (x) * ::Kokkos::real (x) +
                         ::Kokkos::imag (x) * ::Kokkos::imag (x));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return val_type (ArithTraits<mag_type>::zero (), ArithTraits<mag_type>::zero ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return val_type (ArithTraits<mag_type>::one (), ArithTraits<mag_type>::zero ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return val_type (ArithTraits<mag_type>::min (), ArithTraits<mag_type>::min ()); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return val_type (ArithTraits<mag_type>::max (), ArithTraits<mag_type>::max ()); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x.real ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type x) {
    return x.imag ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return ::Kokkos::conj (x);
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const val_type y) {
  //   const mag_type abs_x_square = x.real()*x.real() + x.imag()*x.imag();
  //   const mag_type arg_x = ArithTraits<mag_type>::atan(x.imag()/x.real());
  //   const mag_type half = mag_type(0.5);
  //   const mag_type alpha = (ArithTraits<mag_type>::pow(abs_x_square, half*y.real()) * 
  //                           ArithTraits<mag_type>::exp(-y.imag()*arg_x));
  //   return val_type(alpha* ArithTraits<mag_type>::cos(y.real()*arg_x + half*y.imag()*ArithTraits<mag_type>::log(abs_x_square)),
  //                   alpha* ArithTraits<mag_type>::sin(y.real()*arg_x + half*y.imag()*ArithTraits<mag_type>::log(abs_x_square)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const mag_type y) {
  //   const mag_type arg_x = ArithTraits<mag_type>::atan(x.imag()/x.real());
  //   const mag_type alpha = ArithTraits<mag_type>::pow(abs(x),y);
  //   return val_type(alpha* ArithTraits<mag_type>::cos(y*arg_x), 
  //                   alpha* ArithTraits<mag_type>::sin(y*arg_x));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    return ::Kokkos::sqrt (x);
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
  //   const mag_type r = ::Kokkos::abs(x);
  //   const mag_type phi = ::atan(x.imag()/x.real())/mag_type(3);
  //   const mag_type re = r* ::cos(phi);
  //   const mag_type im = r* ::sin(phi);
  //   return val_type(re,im);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
  //   const mag_type xx = ::exp(x.real());
  //   const mag_type re = xx* ::cos(x.imag());
  //   const mag_type im = xx* ::sin(x.imag());
  //   return val_type(re,im);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
  //   return val_type(ArithTraits<mag_type>::log(abs(x)), ArithTraits<mag_type>::atan(x.imag()/x.real()));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
  //   return log(x)/ArithTraits<mag_type>::log(mag_type(10));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   const val_type ii = val_type(0, 1);
  //   const val_type xx = exp(-ii*x) - exp(ii*x);
  //   const mag_type half = 0.5;
  //   return val_type(-half*xx.imag(),half*xx.real());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   const val_type xx = exp(x) - exp(-x);
  //   const mag_type half = 0.5;
  //   return val_type(half*xx.real(),half*xx.imag());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   const val_type ii(0, 1);
  //   const val_type e_nix = exp(-ii*x);
  //   const val_type e_pix = exp( ii*x);
  //   return ii*(e_nix - e_pix)/(e_nix + e_pix);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   const val_type xx = exp(x) - exp(-x);
  //   const mag_type half = 0.5;
  //   return val_type(half*xx.real(), half*xx.imag());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   const val_type xx = exp(x) + exp(-x);
  //   const mag_type half = 0.5;
  //   return val_type(half*xx.real(), half*xx.imag());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   const val_type e_2x = exp(2*x);
  //   return (e_2x - 1)/(e_2x + 1);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   const val_type ii(0, 1);
  //   const val_type xx = -ii*log(ii*x + sqrt(val_type(1) - x*x));
  //   return val_type(xx.imag(),-xx.real());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   const val_type ii = val_type(0, 1);
  //   const val_type xx = -ii*log(x + ii*sqrt(val_type(1) - x*x));
  //   return val_type(xx.imag(),-xx.real());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   val_type r_val;
  //   const val_type ii = val_type(0, 1);
  //   if (x == ii) {
  //     r_val = val_type(ArithTraits<mag_type>::nan(),  std::numeric_limits<mag_type>::infinity());
  //   } if (x == -ii) {
  //     r_val = val_type(ArithTraits<mag_type>::nan(), -std::numeric_limits<mag_type>::infinity());
  //   } else {
  //     const val_type ii_x = ii*x;
  //     const mag_type half = 0.5;
  //     const val_type xx = log(val_type(1) - ii_x) - log(val_type(1) + ii_x);
  //     r_val = val_type(-half*xx.imag(), half*xx.real());
  //   }
  //   return r_val;
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // ???
    return val_type (ArithTraits<mag_type>::nan (), ArithTraits<mag_type>::nan ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return ArithTraits<mag_type>::epsilon (); // ???
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef ::Kokkos::complex<ArithTraits<mag_type>::halfPrecision> halfPrecision;
  typedef ::Kokkos::complex<ArithTraits<mag_type>::doublePrecision> doublePrecision;

  static const bool isComplex = true;
  static const bool isOrdinal = false;
  static const bool isComparable = false;
  static const bool hasMachineParameters = ArithTraits<mag_type>::hasMachineParameters;
  static bool isnaninf (const val_type& x) {
    return isNan (x) || isInf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static std::string name () {
    return "Kokkos::complex<float>";
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
  //   return sqrt (x);
  // }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type sfmin () {
    return ArithTraits<mag_type>::sfmin (); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION int base () {
    return ArithTraits<mag_type>::base ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type prec () {
    return ArithTraits<mag_type>::prec (); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION int t () {
    return ArithTraits<mag_type>::t ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rnd () {
    return ArithTraits<mag_type>::rnd ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emin () {
    return ArithTraits<mag_type>::emin ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmin () {
    return ArithTraits<mag_type>::rmin ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emax () {
    return ArithTraits<mag_type>::emax ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmax () {
    return ArithTraits<mag_type>::rmax ();
  }
};


template<>
class ArithTraits< ::Kokkos::complex<double> > {
public:
  typedef ::Kokkos::complex<double> val_type;
  typedef double mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = true;

  static constexpr bool has_infinity = true;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() {
    return val_type (ArithTraits<mag_type>::infinity (), ArithTraits<mag_type>::infinity ());
  }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type x) {
    return ArithTraits<mag_type>::isInf (x.real ()) ||
      ArithTraits<mag_type>::isInf (x.imag ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type x) {
    return ArithTraits<mag_type>::isNan (x.real ()) ||
      ArithTraits<mag_type>::isNan (x.imag ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return ::Kokkos::abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return val_type (ArithTraits<mag_type>::zero (), ArithTraits<mag_type>::zero ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return val_type (ArithTraits<mag_type>::one (), ArithTraits<mag_type>::zero ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return val_type (ArithTraits<mag_type>::min (), ArithTraits<mag_type>::min ()); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return val_type (ArithTraits<mag_type>::max (), ArithTraits<mag_type>::max ()); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x.real ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type x) {
    return x.imag ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return ::Kokkos::conj (x);
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const val_type y) {
  //   const mag_type abs_x_square = x.real()*x.real() + x.imag()*x.imag();
  //   const mag_type arg_x = ArithTraits<mag_type>::atan(x.imag()/x.real());
  //   const mag_type half = mag_type(0.5);
  //   const mag_type alpha = (ArithTraits<mag_type>::pow(abs_x_square, half*y.real()) * 
  //                           ArithTraits<mag_type>::exp(-y.imag()*arg_x));
  //   return val_type(alpha* ArithTraits<mag_type>::cos(y.real()*arg_x + half*y.imag()*ArithTraits<mag_type>::log(abs_x_square)),
  //                   alpha* ArithTraits<mag_type>::sin(y.real()*arg_x + half*y.imag()*ArithTraits<mag_type>::log(abs_x_square)));

  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const mag_type y) {
  //   const mag_type arg_x = ArithTraits<mag_type>::atan(x.imag()/x.real());
  //   const mag_type alpha = ArithTraits<mag_type>::pow(abs(x),y);
  //   return val_type(alpha* ArithTraits<mag_type>::cos(y*arg_x), 
  //                   alpha* ArithTraits<mag_type>::sin(y*arg_x));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
     return ::Kokkos::sqrt (x);
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
  //   const mag_type r = ::Kokkos::abs(x);
  //   const mag_type phi = ::atan(x.imag()/x.real())/mag_type(3);
  //   const mag_type re = r* ::cos(phi);
  //   const mag_type im = r* ::sin(phi);
  //   return val_type(re,im);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
  //   const mag_type xx = ::exp(x.real());
  //   const mag_type re = xx* ::cos(x.imag());
  //   const mag_type im = xx* ::sin(x.imag());
  //   return val_type(re,im);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
  //   return val_type(ArithTraits<mag_type>::log(abs(x)), ArithTraits<mag_type>::atan(x.imag()/x.real()));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
  //   return log(x)/ArithTraits<mag_type>::log(mag_type(10));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   const val_type ii = val_type(0, 1);
  //   const val_type xx = exp(-ii*x) - exp(ii*x);
  //   const mag_type half = 0.5;
  //   return val_type(-half*xx.imag(),half*xx.real());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   const val_type xx = exp(x) - exp(-x);
  //   const mag_type half = 0.5;
  //   return val_type(half*xx.real(),half*xx.imag());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   const val_type ii = val_type(0, 1);
  //   const val_type e_nix = exp(-ii*x);
  //   const val_type e_pix = exp( ii*x);
  //   return ii*(e_nix - e_pix)/(e_nix + e_pix);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   const val_type xx = exp(x) - exp(-x);
  //   const mag_type half = 0.5;
  //   return val_type(half*xx.real(), half*xx.imag());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   const val_type xx = exp(x) + exp(-x);
  //   const mag_type half = 0.5;
  //   return val_type(half*xx.real(), half*xx.imag());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   const val_type e_2x = exp(2*x);
  //   return (e_2x - 1)/(e_2x + 1);
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   const val_type ii(0, 1);
  //   const val_type xx = -ii*log(ii*x + sqrt(val_type(1) - x*x));
  //   return val_type(xx.imag(),-xx.real());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   const val_type ii = val_type(0, 1);
  //   const val_type xx = -ii*log(x + ii*sqrt(val_type(1) - x*x));
  //   return val_type(xx.imag(),-xx.real());
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   val_type r_val;
  //   const val_type ii = val_type(0, 1);
  //   if (x == ii) {
  //     r_val = val_type(ArithTraits<mag_type>::nan(),  std::numeric_limits<mag_type>::infinity());
  //   } if (x == -ii) {
  //     r_val = val_type(ArithTraits<mag_type>::nan(), -std::numeric_limits<mag_type>::infinity());
  //   } else {
  //     const val_type ii_x = ii*x;
  //     const mag_type half = 0.5;
  //     const val_type xx = log(val_type(1) - ii_x) - log(val_type(1) + ii_x);
  //     r_val = val_type(-half*xx.imag(), half*xx.real());
  //   }
  //   return r_val;
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // ???
    return val_type (ArithTraits<mag_type>::nan (), ArithTraits<mag_type>::nan ());
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return ArithTraits<mag_type>::epsilon (); // ???
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef ::Kokkos::complex<ArithTraits<mag_type>::halfPrecision> halfPrecision;
  typedef ::Kokkos::complex<ArithTraits<mag_type>::doublePrecision> doublePrecision;

  static const bool isComplex = true;
  static const bool isOrdinal = false;
  static const bool isComparable = false;
  static const bool hasMachineParameters = ArithTraits<mag_type>::hasMachineParameters;
  static bool isnaninf (const val_type& x) {
    return isNan (x) || isInf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static std::string name () {
    return "Kokkos::complex<double>";
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
  //   return sqrt (x);
  // }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type sfmin () {
    return ArithTraits<mag_type>::sfmin (); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION int base () {
    return ArithTraits<mag_type>::base ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type prec () {
    return ArithTraits<mag_type>::prec (); // ???
  }
  static KOKKOS_FORCEINLINE_FUNCTION int t () {
    return ArithTraits<mag_type>::t ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rnd () {
    return ArithTraits<mag_type>::rnd ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emin () {
    return ArithTraits<mag_type>::emin ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmin () {
    return ArithTraits<mag_type>::rmin ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emax () {
    return ArithTraits<mag_type>::emax ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmax () {
    return ArithTraits<mag_type>::rmax ();
  }
};


template<>
class ArithTraits<char> {
public:
  typedef char val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  // The C(++) standard does not require that char be signed.  In
  // fact, signed char, unsigned char, and char are distinct types.
  // We can use std::numeric_limits here because it's a const bool,
  // not a class method.
  static const bool is_signed = std::numeric_limits<char>::is_signed;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    // This avoids warnings based on whether char is signed or unsigned 
    return integer_abs<char>::abs(x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return CHAR_MIN;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return CHAR_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    if (is_signed) {
      return intPowSigned<val_type> (x, y);
    } else {
      return intPowUnsigned<val_type> (x, y);
    }
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    // C++11 defines std::sqrt for integer arguments.  However, we
    // currently can't assume C++11.
    //
    // This cast will result in no loss of accuracy, though it might
    // be more expensive than it should, if we were clever about using
    // bit operations.
    //
    // We take the absolute value first to avoid negative arguments.
    // Negative real arguments to sqrt(float) return (float) NaN, but
    // built-in integer types do not have an equivalent to NaN.
    // Casting NaN to an integer type will thus result in some integer
    // value which appears valid, but is not.  We cannot raise an
    // exception in device functions.  Thus, we prefer to take the
    // absolute value of x first, to avoid issues.  Another
    // possibility would be to test for a NaN output and convert it to
    // some reasonable value (like 0), though this might be more
    // expensive than the absolute value interpreted using the ternary
    // operator.
    return static_cast<val_type> ( ::sqrt (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return static_cast<val_type> ( ::cbrt (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<float> (abs (x))));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<float> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "char";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<signed char> {
public:
  typedef signed char val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x >= 0 ? x : -x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return SCHAR_MIN;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return SCHAR_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowSigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    return static_cast<val_type> ( ::sqrt (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return static_cast<val_type> ( ::cbrt (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<float> (abs (x))));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<float> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "signed char";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<unsigned char> {
public:
  typedef unsigned char val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x; // it's unsigned, so it's positive
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return UCHAR_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowUnsigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<val_type> ( ::sqrt (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return static_cast<val_type> ( ::cbrt (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<float> (x)));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<float> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "unsigned char";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<short> {
public:
  typedef short val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    // std::abs appears to work with CUDA 5.5 at least, but I'll use
    // the ternary expression for maximum generality.  Note that this
    // expression does not necessarily obey the rules for fabs() with
    // NaN input, so it should not be used for floating-point types.
    // It's perfectly fine for signed integer types, though.
    return x >= 0 ? x : -x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    // Macros like this work with CUDA, but
    // std::numeric_limits<val_type>::min() does not, because it is
    // not marked as a __device__ function.
    return SHRT_MIN;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return SHRT_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const val_type y) {
    return intPowSigned<val_type> (x, y);
  }
  //! Integer square root returns a lower bound.
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<val_type> ( ::sqrt (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return static_cast<val_type> ( ::cbrt (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<float> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<float> (abs (x))));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<float> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // short doesn't implement a NaN value, but we can still have it
    // return some "flag" value that can help users find use of
    // uninitialized data.
    return static_cast<val_type> (-1);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "short";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<unsigned short> {
public:
  typedef unsigned short val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x; // it's unsigned, so it's positive
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return USHRT_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowUnsigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<val_type> ( ::sqrt (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return static_cast<val_type> ( ::cbrt (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<float> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<float> (x)));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<float> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<float> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // unsigned short doesn't implement a NaN value, but we can still
    // have it return some "flag" value that can help users find use
    // of uninitialized data.
    return max ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "unsigned short";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<int> {
public:
  typedef int val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    // std::abs appears to work with CUDA 5.5 at least, but I'll use
    // the ternary expression for maximum generality.  Note that this
    // expression does not necessarily obey the rules for fabs() with
    // NaN input, so it should not be used for floating-point types.
    // It's perfectly fine for signed integer types, though.
    return x >= 0 ? x : -x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    // Macros like INT_MIN work with CUDA, but
    // std::numeric_limits<val_type>::min() does not, because it is
    // not marked as a __device__ function.
    return INT_MIN;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return INT_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowSigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<val_type> ( ::sqrt (static_cast<double> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return static_cast<val_type> ( ::cbrt (static_cast<double> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<double> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<double> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<double> (abs (x))));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<double> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // int doesn't implement a NaN value, but we can still have it
    // return some "flag" value that can help users find use of
    // uninitialized data.
    return -1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "int";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<unsigned int> {
public:
  typedef unsigned int val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x; // it's unsigned, so it's positive
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return UINT_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowUnsigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    // This will result in no loss of accuracy, though it might be
    // more expensive than it should, if we were clever about using
    // bit operations.
    return static_cast<val_type> ( ::sqrt (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
    return static_cast<val_type> ( ::cbrt (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<double> (x)));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<double> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // unsigned int doesn't implement a NaN value, but we can still
    // have it return some "flag" value that can help users find use
    // of uninitialized data.
    return max ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "unsigned int";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<long> {
public:
  typedef long val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x >= 0 ? x : -x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return LONG_MIN;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return LONG_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowSigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    using std::sqrt;
    using std::abs;
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    return static_cast<val_type> ( sqrt (static_cast<long double> (abs (x))));
#else
    return static_cast<val_type> ( sqrt (static_cast<double> (abs (x))));
#endif
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<double> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<double> (abs (x))));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<double> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // long doesn't implement a NaN value, but we can still have it
    // return some "flag" value that can help users find use of
    // uninitialized data.
    return -1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "long";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<unsigned long> {
public:
  typedef unsigned long val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return ULONG_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type x, const val_type y) {
    return intPowUnsigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
    using std::sqrt;
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    return static_cast<val_type> ( sqrt (static_cast<long double> (x)));
#else
    return static_cast<val_type> ( sqrt (static_cast<double> (x)));
#endif
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::cbrtl;
    return static_cast<val_type> ( ::cbrtl (static_cast<long double> (x)));
#else
    return static_cast<val_type> ( ::cbrt (static_cast<double> (x)));
#endif
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<long> ( ::log (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<long> ( ::log10 (static_cast<double> (x)));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<double> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // unsigned long doesn't implement a NaN value, but we can still
    // have it return some "flag" value that can help users find use
    // of uninitialized data.
    return max ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "unsigned long";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<long long> {
public:
  typedef long long val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x >= 0 ? x : -x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return LLONG_MIN;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return LLONG_MAX;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowSigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::sqrt;
    using std::abs;
    // IEEE 754 promises that long double has at least 64 significand
    // bits, so we can use it to represent any signed or unsigned
    // 64-bit integer type exactly.  However, CUDA does not implement
    // long double for device functions.
    return static_cast<val_type> ( sqrt (static_cast<long double> (abs (x))));
#else
    // Casting from a 64-bit integer type to double does result in a
    // loss of accuracy.  However, it gives us a good first
    // approximation.  For very large numbers, we may lose some
    // significand bits, but will always get within a factor of two
    // (assuming correct rounding) of the exact double-precision
    // number.  We could then binary search between half the result
    // and twice the result (assuming the latter is <= INT64_MAX,
    // which it has to be, so we don't have to check) to ensure
    // correctness.  It actually should suffice to check numbers
    // within 1 of the result.
    return static_cast<val_type> ( ::sqrt (static_cast<double> (abs (x))));
#endif
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::cbrtl;
    using std::abs;
    return static_cast<val_type> ( cbrtl (static_cast<long double> (abs (x))));
#else
    return static_cast<val_type> ( ::cbrt (static_cast<double> (abs (x))));
#endif
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<double> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<double> (abs (x))));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<double> (abs (x))));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<double> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // long long doesn't implement a NaN value, but we can still have
    // it return some "flag" value that can help users find use of
    // uninitialized data.
    return -1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "long long";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
};


template<>
class ArithTraits<unsigned long long> {
public:
  typedef unsigned long long val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = false;
  static const bool is_integer = true;
  static const bool is_exact = true;
  static const bool is_complex = false;

  static constexpr bool has_infinity = false;
  static KOKKOS_FORCEINLINE_FUNCTION val_type infinity() { return 0; }

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type ) {
    return false;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type x) {
    return x; // unsigned integers are always nonnegative
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return 1;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return ULLONG_MAX ;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag (const val_type ) {
    return 0;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type x) {
    return x;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type
  pow (const val_type x, const val_type y) {
    return intPowUnsigned<val_type> (x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::sqrt;
    return static_cast<val_type> ( sqrt (static_cast<long double> (x)));
#else
    return static_cast<val_type> ( ::sqrt (static_cast<double> (x)));
#endif
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt (const val_type x) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    using std::cbrtl;
    return static_cast<val_type> ( cbrtl (static_cast<long double> (x)));
#else
    return static_cast<val_type> ( ::cbrt (static_cast<double> (x)));
#endif
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type exp (const val_type x) {
    return static_cast<val_type> ( ::exp (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type x) {
    return static_cast<val_type> ( ::log (static_cast<double> (x)));
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type x) {
    return static_cast<val_type> ( ::log10 (static_cast<double> (x)));
  }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
  //   return static_cast<val_type> ( ::sin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
  //   return static_cast<val_type> ( ::cos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
  //   return static_cast<val_type> ( ::tan (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
  //   return static_cast<val_type> ( ::sinh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
  //   return static_cast<val_type> ( ::cosh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
  //   return static_cast<val_type> ( ::tanh (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
  //   return static_cast<val_type> ( ::asin (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
  //   return static_cast<val_type> ( ::acos (static_cast<double> (x)));
  // }
  // static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
  //   return static_cast<val_type> ( ::atan (static_cast<double> (x)));
  // }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    // unsigned long long doesn't implement a NaN value, but we can
    // still have it return some "flag" value that can help users find
    // use of uninitialized data.
    return max ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return zero ();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef val_type halfPrecision;
  typedef val_type doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = true;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude (const val_type x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type x) {
    return conj (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf (const val_type) {
    return false;
  }
  static std::string name () {
    return "unsigned long long";
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type x) {
    return sqrt (x);
  }
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
template<>
struct ArithTraits<dd_real>
{
  typedef dd_real val_type;
  typedef dd_real mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static inline bool isInf (const val_type& x) {
    return isinf (x);
  }
  static inline bool isNan (const val_type& x) {
    return isnan (x);
  }
  static inline mag_type abs (const val_type& x) {
    return ::abs (x);
  }
  static inline val_type zero () {
    return val_type (0.0);
  }
  static inline val_type one () {
    return val_type (1.0);
  }
  static inline val_type min () {
    return std::numeric_limits<val_type>::min ();
  }
  static inline val_type max () {
    return std::numeric_limits<val_type>::max ();
  }
  static inline mag_type real (const val_type& x) {
    return x;
  }
  static inline mag_type imag (const val_type&) {
    return zero ();
  }
  static inline val_type conj (const val_type& x) {
    return x;
  }
  static inline val_type pow (const val_type& x, const val_type& y) {
    return ::pow(x,y);
  }
  static inline val_type sqrt (const val_type& x) {
      return ::sqrt (x);
  }
  static inline val_type cbrt (const val_type& x) {
      return ::cbrt (x);
  }
  static inline val_type exp (const val_type& x) {
      return ::exp (x);
  }
  static inline val_type log (const val_type& x) {
    // dd_real puts its transcendental functions in the global namespace.
    return ::log (x);
  }
  static inline val_type log10 (const val_type& x) {
    return ::log10 (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
    return ::sin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
    return ::cos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
    return ::tan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
    return ::sinh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
    return ::cosh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
    return ::tanh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
    return ::asin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
    return ::acos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
    return ::atan (x);
  }
  static inline val_type nan () {
    return val_type::_nan;
  }
  static val_type epsilon () {
    return std::numeric_limits<val_type>::epsilon ();
  }

  typedef dd_real magnitudeType;
  typedef double halfPrecision;
  typedef qd_real doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;

  static mag_type eps () {
    return epsilon ();
  }
  static mag_type sfmin () {
    return min ();
  }
  static int base ()  {
    return std::numeric_limits<val_type>::radix;
  }
  static mag_type prec () {
    return eps () * base ();
  }
  static int t () {
    return std::numeric_limits<val_type>::digits;
  }
  static mag_type rnd () {
    return std::numeric_limits<val_type>::round_style == std::round_to_nearest ?
      one () :
      zero ();
  }
  static int emin () {
    return std::numeric_limits<val_type>::min_exponent;
  }
  static mag_type rmin () {
    return std::numeric_limits<val_type>::min ();
  }
  static int emax () {
    return std::numeric_limits<val_type>::max_exponent;
  }
  static mag_type rmax () {
    return std::numeric_limits<val_type>::max ();
  }
  static mag_type magnitude (const val_type& x) {
    return ::abs (x);
  }
  static val_type conjugate (const val_type& x) {
    return conj (x);
  }
  static bool isnaninf (const val_type& x) {
    return isNan (x) || isInf (x);
  }
  static std::string name () { return "dd_real"; }
  static val_type squareroot (const val_type& x) {
    return ::sqrt (x);
  }
};


template<>
struct ArithTraits<qd_real>
{
  typedef qd_real val_type;
  typedef qd_real mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static inline bool isInf (const val_type& x) {
    return isinf (x);
  }
  static inline bool isNan (const val_type& x) {
    return isnan (x);
  }
  static inline mag_type abs (const val_type& x) {
    return ::abs (x);
  }
  static inline val_type zero () {
    return val_type (0.0);
  }
  static inline val_type one () {
    return val_type (1.0);
  }
  static inline val_type min () {
    return std::numeric_limits<val_type>::min ();
  }
  static inline val_type max () {
    return std::numeric_limits<val_type>::max ();
  }
  static inline mag_type real (const val_type& x) {
    return x;
  }
  static inline mag_type imag (const val_type&) {
    return zero ();
  }
  static inline val_type conj (const val_type& x) {
    return x;
  }
  static inline val_type pow (const val_type& x, const val_type& y) {
    return ::pow (x, y);
  }
  static inline val_type sqrt (const val_type& x) {
    return ::sqrt (x);
  }
  static inline val_type cbrt (const val_type& x) {
    return ::cbrt (x);
  }
  static inline val_type exp (const val_type& x) {
    return ::exp (x);
  }
  static inline val_type log (const val_type& x) {
    // val_type puts its transcendental functions in the global namespace.
    return ::log (x);
  }
  static inline val_type log10 (const val_type& x) {
    return ::log10 (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sin (const val_type x) {
    return ::sin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cos (const val_type x) {
    return ::cos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tan (const val_type x) {
    return ::tan (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sinh (const val_type x) {
    return ::sinh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type cosh (const val_type x) {
    return ::cosh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type tanh (const val_type x) {
    return ::tanh (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type asin (const val_type x) {
    return ::asin (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type acos (const val_type x) {
    return ::acos (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type atan (const val_type x) {
    return ::atan (x);
  }
  static inline val_type nan () {
    return val_type::_nan;
  }
  static inline val_type epsilon () {
    return std::numeric_limits<val_type>::epsilon ();
  }

  typedef qd_real magnitudeType;
  typedef dd_real halfPrecision;
  // The QD library does not have an "oct-double real" class.  One
  // could use an arbitrary-precision library like MPFR or ARPREC,
  // with the precision set appropriately, to get an
  // extended-precision type for qd_real.
  typedef qd_real doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;

  static mag_type eps () {
    return epsilon ();
  }
  static mag_type sfmin () {
    return min ();
  }
  static int base ()  {
    return std::numeric_limits<val_type>::radix;
  }
  static mag_type prec () {
    return eps () * base ();
  }
  static int t () {
    return std::numeric_limits<val_type>::digits;
  }
  static mag_type rnd () {
    return std::numeric_limits<val_type>::round_style == std::round_to_nearest ?
      one () :
      zero ();
  }
  static int emin () {
    return std::numeric_limits<val_type>::min_exponent;
  }
  static mag_type rmin () {
    return std::numeric_limits<val_type>::min ();
  }
  static int emax () {
    return std::numeric_limits<val_type>::max_exponent;
  }
  static mag_type rmax () {
    return std::numeric_limits<val_type>::max ();
  }
  static mag_type magnitude (const val_type& x) {
    return ::abs (x);
  }
  static val_type conjugate (const val_type& x) {
    return conj (x);
  }
  static bool isnaninf (const val_type& x) {
    return isNan (x) || isInf (x);
  }
  static std::string name () { return "qd_real"; }
  static val_type squareroot (const val_type& x) {
    return ::sqrt (x);
  }
};
#endif // HAVE_KOKKOS_QD

} // namespace Details

  // Promote ArithTraits into Kokkos namespace.  At some point, we
  // will remove it from the Details namespace completely.  We leave
  // it there for now, because a lot of code depends on it being
  // there.
  using Details::ArithTraits;
} // namespace Kokkos

#endif // KOKKOS_ARITHTRAITS_HPP
