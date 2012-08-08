// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_AS_HPP
#define TEUCHOS_AS_HPP

#include "Teuchos_Assert.hpp"
#include <limits>

#ifdef HAVE_TEUCHOS_QD
#include <qd/qd_real.h>
#include <qd/dd_real.h>
#endif

namespace Teuchos {


/** \class ValueTypeConversionTraits
 * \brief Default traits class for all conversions between value types.
 * \ingroup teuchos_language_support_grp
 *
 * \note Users should never call this class directly.  Please use the
 *   <tt>as()</tt> or <tt>asSafe()</tt> template functions.
 *
 * \tparam TypeTo The type to which to convert; the output type.
 * \tparam TypeFrom The type from which to convert; the input type.
 *
 * The type conversion functions <tt>as()</tt> and <tt>asSafe()</tt>
 * use this traits class to convert between types.  The default
 * implementation of this class simply does an implicit type
 * conversion.  Syntactically, it expects either that TypeTo have a
 * constructor which takes a single TypeFrom argument, like this:
 * \code
 * class TypeTo {
 * public:
 *   TypeTo (const TypeFrom& x);
 * }
 * \endcode
 * or that TypeFrom has an "operator TypeTo()" method, like this:
 * \code
 * class TypeFrom {
 * public:
 *   operator TypeTo ();
 * }
 * \endcode
 * Any conversions which are built into C++ and are safe to do will
 * not need a traits class specialization, and should not generate any
 * compiler warnings.  This includes the conversions <tt>float</tt> to
 * <tt>double</tt>, <tt>short type</tt> to <tt>type</tt>,
 * <tt>type</tt> to <tt>long type</tt>, or an enum value to
 * <tt>int</tt> (where <tt>type</tt> may be either <tt>int</tt> or
 * <tt>unsigned int</tt>).
 *
 * Any conversion which is not syntactically legal according to the
 * above rules _must_ have a specialization.  There are a number of
 * examples in the header file below, including <tt>qd_real</tt>.
 * Other examples include <tt>std::string</tt> to <tt>int</tt>,
 * <tt>double</tt>, etc.
 *
 * Any conversion that _is_ syntactically legal, but could cause
 * compiler warnings and/or result in incorrect behavior at run time
 * should be given a specialization that does not rely on an implicit
 * conversion.  This includes the following:
 *
 * - <tt>type</tt> to and from <tt>unsigned type</tt>, where
 *   <tt>type</tt> is a built-in integer type
 * - <tt>double</tt> to <tt>int</tt>, or between any floating-point
 *   and integer types where overflow is possible
 *
 * If the user (through <tt>as()</tt> or <tt>asSafe()</tt>) requests a
 * conversion for which no specialization of this class exists), then
 * the default implementation below will be instantiated.  If the
 * conversion is not syntactically correct, then the compiler will
 * report a compiler error.  If the conversion is syntactically
 * correct but unsafe, the compiler _may_ report a warning.  In either
 * case, you can fix the error or warning by specializing this class
 * for your combination of types.  There are a number of examples of
 * specializations in this header file, some of which include bounds
 * checking for overflow (for safeConvert()).
 *
 * \note We cannot promise that converting from T1 to T2 and back
 *   again will result in the same T1 value with which we started.
 *   For example, converting from a long long to a double may result
 *   in truncation, since long long has 63 bits of significand and
 *   double has 53.
 */
template<class TypeTo, class TypeFrom>
class ValueTypeConversionTraits {
public:
  //! Convert t from a TypeFrom object to a TypeTo object.
  static TypeTo convert (const TypeFrom t) {
    // This default implementation is just an implicit conversion and
    // may generate compiler warnings on dangerous conversions.
    return t;
  }

  //! Convert t from a TypeFrom object to a TypeTo object, in a more safe way.
  static TypeTo safeConvert (const TypeFrom t) {
    // This default implementation is just an implicit conversion and
    // may generate compiler warnings on dangerous conversions.  No
    // runtime checking (e.g., for overflow) can be done by default;
    // only specializations can define meaningful and portable
    // run-time checks of conversions.
    return t;
  }
};

/** \fn as
 * \brief Convert from one value type to another.
 * \ingroup teuchos_language_support_grp
 *
 * User documentation
 * ==================
 *
 * This template function lets you convert from one value type to
 * another, possibly with checks for overflow (where appropriate) in a
 * debug build of Teuchos.  For example, to convert between int and
 * double:
 * \code
 * double d = 3.14;
 * int i = Teuchos::as<double> (d);
 * assert (i == 3);
 * \endcode
 * In a debug build of Teuchos, this will check for overflow, since
 * there are some double-precision floating-point values that do not
 * fit in a 32-bit integer.  In a release build, this will not check
 * for overflow.  You are responsible for knowing the difference.  If
 * you _always_ want to check for overflow (e.g., to validate user
 * input), use the asSafe() function.  Note that conversion from a
 * floating-point number to an integer, or from a higher-precision
 * floating-point number to a lower-precision floating-point number,
 * may truncate or round (as it does in the above example).
 *
 * "Debug build of Teuchos" means more than just building with debug
 * compiler flags.  It means debug checking was turned on when
 * Trilinos was built.  If you are building Trilinos yourself, you may
 * turn on debug checking by setting the Trilinos_ENABLE_DEBUG CMake
 * configure option to ON (rather than OFF, which is the default).
 * Note that enabling debug checking affects other operations in
 * Teuchos besides this conversion, and may have a significant
 * run-time cost, especially for RCP and ArrayRCP.
 *
 * \note We cannot promise that converting from a type T1 to another
 *   type T2 and back again will result in the same T1 value with
 *   which we started.  For example, converting from a long long to a
 *   double may result in truncation, since long long has 63 bits of
 *   significand and double has 53.
 *
 * Developer documentation
 * =======================
 *
 * This function just uses the traits class ValueTypeConversionTraits
 * to perform the actual conversion.  If debug checking is turned on,
 * this function uses the traits class' safeConvert() method to
 * perform possibly checked conversion.  Otherwise, it uses the traits
 * class' convert() method for unchecked conversion.
 *
 * If you want to specialize this function's behavior, you should
 * specialize ValueTypeConversionTraits for your combination of input
 * and output types (TypeFrom resp. TypeTo).  Be sure to define the
 * specialization in the Teuchos namespace.  We provide
 * specializations of ValueTypeConversionTraits for a variety of
 * types.  You must define both safeConvert() and convert() in the
 * specialization, since as() will call safeConvert() in a debug build
 * and convert() in a release build.
 *
 * \note The implementations below do not consider truncation of
 *   floating-point values to be unsafe conversion.  For example,
 *   converting from a long long (63 bits of significand) to a double
 *   (53 bits of significand) may result in truncation, but we do not
 *   consider this unsafe.  "Unsafe" mainly refers to overflow or lack
 *   of representation.
 */
template<class TypeTo, class TypeFrom>
inline TypeTo as( const TypeFrom& t )
{
#ifdef HAVE_TEUCHOS_DEBUG
  return ValueTypeConversionTraits<TypeTo,TypeFrom>::safeConvert(t);
#else
  return ValueTypeConversionTraits<TypeTo,TypeFrom>::convert(t);
#endif // HAVE_TEUCHOS_DEBUG
}


/** \fn asSafe
 * \brief Convert from one value type to another,
 *   checking for validity first if appropriate.
 * \ingroup teuchos_language_support_grp
 *
 * User documentation
 * ==================
 *
 * This template function lets you convert from one value type to
 * another.  For example, to convert between int and double:
 * \code
 * double d = 3.14;
 * int i = Teuchos::as<double> (d);
 * assert (i == 3);
 * \endcode
 * This function always checks for validity of the conversion before
 * attempting it.  Examples of invalid conversions are those which
 * would overflow, for example from double to int, if the
 * double-precision floating-point input is bigger than the largest
 * int or smaller than the smallest int.  If you only which to check
 * in a debug build, use the as() function instead.  Note that
 * conversion from a floating-point number to an integer, or from a
 * higher-precision floating-point number to a lower-precision
 * floating-point number, may truncate or round (as it does in the
 * above example).
 *
 * Developer documentation
 * =======================
 *
 * This function just uses the traits class ValueTypeConversionTraits
 * to perform the actual conversion.  It always uses the traits class'
 * safeConvert() method to perform a possibly checked conversion.
 *
 * If you want to specialize this function's behavior, you should
 * specialize ValueTypeConversionTraits for your combination of input
 * and output types (TypeFrom resp. TypeTo).  Be sure to define the
 * specialization in the Teuchos namespace.  We provide
 * specializations of ValueTypeConversionTraits for a variety of
 * types.
 *
 * \note The implementations below do not consider truncation of
 *   floating-point values to be unsafe conversion.  For example,
 *   converting from a long long (63 bits of significand) to a double
 *   (53 bits of significand) may result in truncation, but we do not
 *   consider this unsafe.  "Unsafe" mainly refers to overflow or lack
 *   of representation.
 */
template<class TypeTo, class TypeFrom>
inline TypeTo asSafe( const TypeFrom& t )
{
  return ValueTypeConversionTraits<TypeTo,TypeFrom>::safeConvert(t);
}


/// \class asFunc
/// \brief Function object wrapper for as().
/// \ingroup teuchos_language_support_grp
///
/// \tparam TypeTo Type to which to convert; the output type.
///
/// Sometimes it is useful to pass around the as() type conversion
/// function as a first-class object, for example as a function
/// argument of generic algorithms such as std::transform().  In this
/// case, you may use this class, which invokes as() in its operator()
/// method.
template <class TypeTo>
class asFunc {
  public:
  asFunc() {}

  template <class TypeFrom>
  inline TypeTo operator()(const TypeFrom &t) {
    return as<TypeTo>(t);
  }
};


//
// Standard specializations of ValueTypeConversionTraits
//

//! Convert double to int, without compiler warnings, with optional range check.
template<>
class ValueTypeConversionTraits<int, double> {
public:
  /// \brief Convert the given double to an int.
  ///
  /// \warning Double-precision floating-point values (64 bits) may
  ///   overflow int (32 bits).  You should use safeConvert()
  ///   if you aren't sure that the given value fits in int.
  static int convert (const double t) {
    // Implicit conversion from double to int causes compiler
    // warnings, but static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert the given double to an int, checking for overflow first.
  static int safeConvert (const double t) {
    const int minInt = std::numeric_limits<int>::min ();
    const int maxInt = std::numeric_limits<int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<double> (minInt) || t > static_cast<double> (maxInt),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<int, double>::safeConvert: "
      "Input double t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to int.");
    return static_cast<int> (t);
  }
};

//! Convert float to int, without compiler warnings.
template<>
class ValueTypeConversionTraits<int, float> {
public:
  //! Convert the given float to an int.
  static int convert (const float t) {
    // Implicit conversion from float to int may cause compiler
    // warnings, but static_cast does not.  Note that casting from
    // (32-bit) float to (32-bit signed) int can never overflow.
    return static_cast<int> (t);
  }

  //! Convert the given float to an int.
  static int safeConvert (const float t) {
    // Implicit conversion from float to int may cause compiler
    // warnings, but static_cast does not.  Note that casting from
    // (32-bit) float to (32-bit signed) int can never overflow.
    return static_cast<int> (t);
  }
};


#ifdef HAVE_TEUCHOS_LONG_LONG_INT

//! Convert long long to float, without compiler warnings, with optional overflow check.
template<>
class ValueTypeConversionTraits<float, long long> {
public:
  /// \brief Convert the given long long to a float.
  ///
  /// \warning long long (64 bits) integer values may overflow float
  ///   (32 bits).  You should use safeConvert() if you aren't sure
  ///   that the given value fits a float.
  static float convert (const long long t) {
    // Implicit conversion from long long to float may cause compiler
    // warnings, but static_cast does not.
    return static_cast<float> (t);
  }

  //! Convert from long long to float, checking for overflow first.
  static float safeConvert (const long long t) {
    const float minFloat = std::numeric_limits<float>::min ();
    const float maxFloat = std::numeric_limits<float>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<long long> (minFloat) ||
      t > static_cast<long long> (maxFloat),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<float, long long>::safeConvert: "
      "Input long long t = " << t << " is out of the valid range [" << minFloat
      << ", " << maxFloat << "] for conversion to float.");

    // Implicit conversion from long long to float may cause compiler
    // warnings, but static_cast does not.
    return static_cast<float> (t);
  }
};

//! Convert unsigned long long to float, without compiler warnings, with optional overflow check.
template<>
class ValueTypeConversionTraits<float, unsigned long long> {
public:
  /// \brief Convert the given unsigned long long to a float.
  ///
  /// \warning unsigned long long (64 bits) integer values may
  ///   overflow float (32 bits).  You should use safeConvert() if you
  ///   aren't sure that the given value fits a float.
  static float convert (const unsigned long long t) {
    // Implicit conversion from unsigned long long to float may cause
    // compiler warnings, but static_cast does not.
    return static_cast<float> (t);
  }

  //! Convert from unsigned long long to float, checking for overflow first.
  static float safeConvert (const unsigned long long t) {
    const float minFloat = std::numeric_limits<float>::min ();
    const float maxFloat = std::numeric_limits<float>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<unsigned long long> (minFloat) ||
      t > static_cast<unsigned long long> (maxFloat),
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<float, unsigned long long>::safeConvert: "
      "Input unsigned long long t = " << t << " is out of the valid range [" << minFloat
      << ", " << maxFloat << "] for conversion to float.");

    // Implicit conversion from unsigned long long to float may cause
    // compiler warnings, but static_cast does not.
    return static_cast<float> (t);
  }
};

//! Convert long long to int, without compiler warnings, with optional overflow check.
template<>
class ValueTypeConversionTraits<int, long long> {
public:
  /// \brief Convert the given long long to a int.
  ///
  /// \warning long long (64 bits) integer values may overflow int
  ///   (32 bits).  You should use safeConvert() if you aren't sure
  ///   that the given value fits a int.
  static int convert (const long long t) {
    // Implicit conversion from long long to int may cause compiler
    // warnings, but static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert from long long to int, checking for overflow first.
  static int safeConvert (const long long t) {
    const int minInt = std::numeric_limits<int>::min ();
    const int maxInt = std::numeric_limits<int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<long long> (minInt) ||
      t > static_cast<long long> (maxInt),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<int, long long>::safeConvert: "
      "Input long long t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to int.");

    // Implicit conversion from long long to int may cause compiler
    // warnings, but static_cast does not.
    return static_cast<int> (t);
  }
};


//! Convert long long to unsigned int, without compiler warnings, with optional overflow check.
template<>
class ValueTypeConversionTraits<unsigned int, long long> {
public:
  /// \brief Convert the given long long to a unsigned int.
  ///
  /// \warning long long (64 bits) integer values may overflow
  ///   unsigned int (32 bits).  You should use safeConvert() if you
  ///   aren't sure that the given value fits a unsigned int.
  static unsigned int convert (const long long t) {
    // Implicit conversion from long long to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert from long long to unsigned int, checking for overflow first.
  static unsigned int safeConvert (const long long t) {
    const unsigned int minInt = std::numeric_limits<unsigned int>::min ();
    const unsigned int maxInt = std::numeric_limits<unsigned int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<long long> (minInt) ||
      t > static_cast<long long> (maxInt),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<unsigned int, long long>::safeConvert: "
      "Input long long t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to unsigned int.");

    // Implicit conversion from long long to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }
};

//! Convert unsigned long long to int, without compiler warnings, with optional overflow check.
template<>
class ValueTypeConversionTraits<int, unsigned long long> {
public:
  /// \brief Convert the given unsigned long long to a int.
  ///
  /// \warning unsigned long long (64 bits) integer values may
  ///   overflow int (32 bits).  You should use safeConvert() if you
  ///   aren't sure that the given value fits a int.
  static int convert (const unsigned long long t) {
    // Implicit conversion from unsigned long long to int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert from unsigned long long to int, checking for overflow first.
  static int safeConvert (const unsigned long long t) {
    const int minInt = std::numeric_limits<int>::min ();
    const int maxInt = std::numeric_limits<int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<unsigned long long> (minInt) ||
      t > static_cast<unsigned long long> (maxInt),
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<int, unsigned long long>::safeConvert: "
      "Input unsigned long long t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to int.");

    // Implicit conversion from unsigned long long to int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<int> (t);
  }
};

//! Convert unsigned long long to unsigned int, without compiler warnings, with optional overflow check.
template<>
class ValueTypeConversionTraits<unsigned int, unsigned long long> {
public:
  /// \brief Convert the given unsigned long long to a unsigned int.
  ///
  /// \warning unsigned long long (64 bits) integer values may
  ///   overflow unsigned int (32 bits).  You should use safeConvert()
  ///   if you aren't sure that the given value fits a unsigned int.
  static unsigned int convert (const unsigned long long t) {
    // Implicit conversion from unsigned long long to unsigned int may
    // cause compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert from unsigned long long to unsigned int, checking for overflow first.
  static unsigned int safeConvert (const unsigned long long t) {
    const unsigned int minInt = std::numeric_limits<unsigned int>::min ();
    const unsigned int maxInt = std::numeric_limits<unsigned int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<unsigned long long> (minInt) ||
      t > static_cast<unsigned long long> (maxInt),
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<unsigned int, unsigned long long>::safeConvert: "
      "Input unsigned long long t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to unsigned int.");

    // Implicit conversion from unsigned long long to unsigned int may
    // cause compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }
};

#endif // HAVE_TEUCHOS_LONG_LONG_INT


/** \brief Convert raw C string to std::string. */
template<int N>
class ValueTypeConversionTraits<std::string, char[N]> {
public:
  static std::string convert( const char t[] )
    { return std::string(t); }
  static std::string safeConvert( const char t[] )
    { return std::string(t); }
};

#ifdef HAVE_TEUCHOS_QD

/** \brief Convert qd_real to double. */
template <>
class ValueTypeConversionTraits<double, qd_real> {
public:
  inline static double convert( const qd_real t )
    { return to_double(t); }
  inline static double safeConvert( const qd_real t )
    { return to_double(t); }
};

/** \brief Convert qd_real to float. */
template <>
class ValueTypeConversionTraits<float, qd_real> {
public:
  inline static float convert( const qd_real t )
    { return (float)to_double(t); }
  inline static float safeConvert( const qd_real t )
    { return (float)to_double(t); }
};

/** \brief Convert qd_real to int. */
template <>
class ValueTypeConversionTraits<int, qd_real> {
public:
  inline static int convert( const qd_real t )
    { return to_int(t); }
  inline static int safeConvert( const qd_real t )
    { return to_int(t); }
};

/** \brief Convert qd_real to dd_real. */
template <>
class ValueTypeConversionTraits<dd_real, qd_real> {
public:
  inline static dd_real convert( const qd_real t )
    { return to_dd_real(t); }
  inline static dd_real safeConvert( const qd_real t )
    { return to_dd_real(t); }
};

/** \brief Convert dd_real to double. */
template <>
class ValueTypeConversionTraits<double, dd_real> {
public:
  inline static double convert( const dd_real t )
    { return to_double(t); }
  inline static double safeConvert( const dd_real t )
    { return to_double(t); }
};

/** \brief Convert dd_real to float. */
template <>
class ValueTypeConversionTraits<float, dd_real> {
public:
  inline static float convert( const dd_real t )
    { return (float)to_double(t); }
  inline static float safeConvert( const dd_real t )
    { return (float)to_double(t); }
};

/** \brief Convert dd_real to int. */
template <>
class ValueTypeConversionTraits<int, dd_real> {
public:
  inline static int convert( const dd_real t )
    { return to_int(t); }
  inline static int safeConvert( const dd_real t )
    { return to_int(t); }
};

/** \brief Convert long unsigned int to qd_real. */
template <>
class ValueTypeConversionTraits<qd_real, long unsigned int> {
public:
  inline static qd_real convert( const long unsigned int t )
    { return ValueTypeConversionTraits<qd_real,int>::convert(ValueTypeConversionTraits<int,long unsigned int>::convert(t)); }
  inline static qd_real safeConvert( const long unsigned int t )
    { return ValueTypeConversionTraits<qd_real,int>::safeConvert(ValueTypeConversionTraits<int,long unsigned int>::safeConvert(t)); }
};

/** \brief Convert long unsigned int to dd_real. */
template <>
class ValueTypeConversionTraits<dd_real, long unsigned int> {
public:
  inline static dd_real convert( const long unsigned int t )
    { return ValueTypeConversionTraits<dd_real,int>::convert(ValueTypeConversionTraits<int,long unsigned int>::convert(t)); }
  inline static dd_real safeConvert( const long unsigned int t )
    { return ValueTypeConversionTraits<dd_real,int>::safeConvert(ValueTypeConversionTraits<int,long unsigned int>::safeConvert(t)); }
};

#endif // HAVE_TEUCHOS_QD

// ToDo: Add more specializations as needed!


} // end namespace Teuchos


#endif // TEUCHOS_AS_HPP
