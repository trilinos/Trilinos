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

/// \file Teuchos_as.hpp
/// \brief Definition of Teuchos::as, for conversions between types.
///
/// This header file declares Teuchos::as, which is a template
/// function that converts between two different types.  For example,
/// the following code converts from \c double to \c int:
/// \code
/// double d = 3.14;
/// int i = Teuchos::as<double> (d);
/// assert (i == 3);
/// \endcode
/// In a debug build of Teuchos, this code would check for overflow
/// when converting from double (64 bits) to int (32 bits, on most
/// platforms these days), and throw an exception if overflow occurs.
/// In a release build, this code would not check for overflow.  Users
/// who want to check for overflow may use the Teuchos::asSafe
/// template function.  This works the same way as Teuchos::as, except
/// that it always checks for overflow.
///
/// This file includes definitions of a small number of useful
/// conversions.  If you want to define your own conversions, you
/// should specialize the Teuchos::ValueTypeConversionTraits template
/// class for the "to" and "from" types in your conversion.  This
/// automatically defines Teuchos::as and Teuchos::asSafe for your
/// conversion.

#include "Teuchos_Assert.hpp"
#include <limits>
#include <cstdlib>
#include <cerrno>
#include <climits>

#ifdef HAVE_TEUCHOS_QD
#include <qd/qd_real.h>
#include <qd/dd_real.h>
#endif // HAVE_TEUCHOS_QD

namespace Teuchos {


/** \class ValueTypeConversionTraits
 * \brief Default traits class for all conversions between value types.
 * \ingroup teuchos_language_support_grp
 *
 * \note Users should never call this class directly.  Please use the
 *   as() or asSafe() template functions.
 *
 * \tparam TypeTo The type to which to convert; the output type.
 * \tparam TypeFrom The type from which to convert; the input type.
 *
 * The type conversion functions as() and asSafe() defined in this
 * file use this traits class to convert between types.  The default
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
 * <tt>double</tt>, <tt>short</tt> to <tt>int</tt>, <tt>int</tt> to
 * <tt>long</tt>, or an enum value to <tt>int</tt>.
 *
 * Any conversion which is not syntactically legal according to the
 * above rules <i>must</i> have a specialization.  There are a number
 * of examples in the header file below, including conversions between
 * <tt>qd_real</tt> and various built-in types for which
 * <tt>qd_real</tt> does not provide a native conversion operator.
 * Other examples include <tt>std::string</tt> to <tt>int</tt> or
 * <tt>double</tt>.
 *
 * Any conversion that <i>is</i> syntactically legal, but could cause
 * compiler warnings and/or result in incorrect behavior at run time
 * should be given a specialization that does not rely on an implicit
 * conversion.  This includes the following:
 *
 * - <tt>type</tt> to and from <tt>unsigned type</tt>, where
 *   <tt>type</tt> is a built-in integer type
 * - <tt>double</tt> to <tt>int</tt>, or between any floating-point
 *   and integer types where overflow is possible
 *
 * If the user (through as() or asSafe()) requests a conversion for
 * which no specialization of this class exists, then the default
 * implementation below will be instantiated.  If the conversion is
 * not syntactically correct, then the compiler will report a compiler
 * error.  If the conversion is syntactically correct but unsafe, the
 * compiler <i>may</i> report a warning.  In either case, you can fix
 * the error or warning by specializing this class for your
 * combination of types.  There are a number of examples of
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

  //! Convert t from a TypeFrom object to a TypeTo object, with checks for validity.
  static TypeTo safeConvert (const TypeFrom t) {
    // This default implementation is just an implicit conversion and
    // may generate compiler warnings on dangerous conversions.  No
    // runtime checking (e.g., for overflow) can be done by default;
    // only specializations can define meaningful and portable
    // run-time checks of conversions.
    return t;
  }
};

/** \brief Convert from one value type to another.
 * \ingroup teuchos_language_support_grp
 *
 * \tparam TypeTo The type to which to convert; the output type.
 * \tparam TypeFrom The type from which to convert; the input type.
 *
 * \section Teuchos_as_User User documentation
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
 * you <i>always</i> want to check for overflow (e.g., to validate
 * user input), use the asSafe() function.  Note that conversion from
 * a floating-point number to an integer, or from a higher-precision
 * floating-point number to a lower-precision floating-point number,
 * may truncate or round (as it does in the above example).  We do
 * not check for truncation or rounding.
 *
 * "Debug build of Teuchos" means more than just building with debug
 * compiler flags.  It means debug checking was turned on when
 * Trilinos was built.  If you are building Trilinos yourself, you may
 * turn on debug checking by setting the
 * <tt>Trilinos_ENABLE_DEBUG</tt> CMake configure option to \c ON
 * (rather than \c OFF, which is the default).  Note that enabling
 * debug checking affects other operations in Teuchos besides this
 * conversion, and may have a significant run-time cost, especially
 * for RCP and ArrayRCP.
 *
 * \note We cannot promise that converting from a type T1 to another
 *   type T2 and back again will result in the same T1 value with
 *   which we started.  For example, converting from a long long to a
 *   double may result in truncation, since long long has 63 bits of
 *   significand and double has 53.
 *
 * \section Teuchos_as_Dev Developer documentation
 *
 * This function uses the traits class ValueTypeConversionTraits to
 * perform checking and conversion.  If debug checking is turned on,
 * this function uses the traits class' safeConvert() method to
 * perform possibly checked conversion.  Otherwise, it uses the traits
 * class' convert() method for unchecked conversion.
 *
 * If you want to specialize this function's behavior, you should
 * specialize ValueTypeConversionTraits for your combination of input
 * and output types (TypeFrom resp. TypeTo).  Be sure to define the
 * specialization in the Teuchos namespace.  We provide
 * specializations of ValueTypeConversionTraits in this file for a
 * variety of types.  You must define both safeConvert() and convert()
 * in the specialization, since as() will call safeConvert() in a
 * debug build and convert() in a release build.
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


/** \brief Convert from one value type to another, with validity checks if appropriate.
 * \ingroup teuchos_language_support_grp
 *
 * \tparam TypeTo The type to which to convert; the output type.
 * \tparam TypeFrom The type from which to convert; the input type.
 *
 * \section Teuchos_asSafe_User User documentation
 *
 * This template function lets you convert from one value type to
 * another.  For example, to convert between int and double:
 * \code
 * double d = 3.14;
 * int i = Teuchos::asSafe<double> (d);
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
 * \section Teuchos_asSafe_Dev Developer documentation
 *
 * This function uses the traits class ValueTypeConversionTraits to
 * perform the actual checking and conversion.  It always uses the
 * traits class' safeConvert() method to perform a possibly checked
 * conversion.
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
/// method.  The operator() method is templated on the input type
/// TypeFrom.
template <class TypeTo>
class asFunc {
public:
  asFunc() {}

  template <class TypeFrom>
  inline TypeTo operator()(const TypeFrom &t) {
    return as<TypeTo>(t);
  }
};



namespace { // anonymous

  /// Helper function for converting the given \c std::string to an
  /// integer of type <tt>IntType</tt>.
  ///
  /// \tparam IntType A built-in integer type, like \c int or \c long.
  ///   It may be signed or unsigned.
  ///
  /// \param t [in] The string to convert.
  ///
  /// \param rawConvert [in] A function with the same arguments as
  ///   strtol, strtoul, strtoll, or strtoull, which returns
  ///   <tt>IntType<tt>.  It must return the same type as
  ///   <tt>IntType</tt>.  Note that all of these but strtol require
  ///   C99 support.  It's up to you to pick the right function for
  ///   <tt>IntType</tt>.
  ///
  /// \param intTypeName [in] Human-readable string which is the name
  ///   of <tt>IntType</tt>.
  template<class IntType>
  IntType
  intToString (const std::string& t,
               IntType (*rawConvert) (const char*, char**, int),
               const char* intTypeName)
  {
    // We call the "raw" conversion function instead of using
    // std::istringstream, because we want more detailed information
    // in case of failure to convert.  I have no idea what
    // operator>>(std::istream&, unsigned long long&) does if it
    // encounters an integer too long to fit in IntType, for example.
    //
    // mfh 13 Nov 2012: It's fair to assume that if you have "long
    // long", then your implementation of the C standard library
    // includes strtoul().  Ditto for "unsigned long long" and
    // strtoull().  If this is not the case, we could include a
    // configure-time test for these functions(), with a fall-back to
    // an std::istringstream operator>> implementation.
    char* endptr = NULL;
    // Keep the pointer, because std::string doesn't necessarily
    // guarantee that this is the same across calls to c_str(), does
    // it?  Or perhaps it does...
    const char* t_ptr = t.c_str ();
    // We preset errno to 0, to distinguish success or failure after
    // calling strtoull.  Most implementations of the C standard
    // library written with threads in mind have errno be a macro that
    // expands to thread-local storage.  Thanks to the Linux
    // documentation for strtol ("man 3 strtol", Red Hat Enterprise
    // Linux 5) for advice with the following checks.
    errno = 0;
    const IntType val = rawConvert (t_ptr, &endptr, 10);

    const IntType minVal = std::numeric_limits<IntType>::min ();
    const IntType maxVal = std::numeric_limits<IntType>::max ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      errno == ERANGE && (val == minVal || val == maxVal),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<" << intTypeName << ", std::string>::convert: "
      "The integer value in the given string \"" << t << "\" overflows " << intTypeName << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      errno != 0 && val == 0,
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<" << intTypeName << ", std::string>::convert: "
      "The conversion function was unable to convert the given string \"" << t << "\" to " << intTypeName << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      endptr == t_ptr, // See above discussion of c_str().
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<" << intTypeName << ", std::string>::convert: "
      "The conversion function was unable to read any integer digits from the given string "
      "\"" << t << "\".");
    return val;
  }

  /// Helper function for converting the given \c std::string to a
  /// real-valued floating-point number of type RealType.
  ///
  /// \tparam RealType The real-valued floating-point type of the
  ///   return value.  We always assume that RealType is default
  ///   constructible and that <tt>operator>>(std::istream&,
  ///   RealType&)</tt> is defined.  We also assume that <tt>t !=
  ///   0</tt> is a well-formed Boolean expression for t of type
  ///   RealType.
  ///
  /// \param t [in] The string to convert.
  ///
  /// \param rawConvert [in] If not NULL, this must be one of the
  ///   following C standard library functions: strtod, strtof, or
  ///   strtold.  (strtof and strtold require C99 support.)  In that
  ///   case, we use this function to read the value from the input
  ///   string.  If NULL, we use <tt>operator>>(std::istream&,
  ///   RealType&)</tt> to read the value from the string (via
  ///   std::istringstream).
  ///
  /// \param realTypeName [in] Human-readable string which is the name
  ///   of RealType.
  template<class RealType>
  RealType
  realToString (const std::string& t,
               RealType (*rawConvert) (const char*, char**),
               const char* realTypeName)
  {
    if (rawConvert == NULL) {
      std::istringstream in (t);
      RealType out;
      in >> out;
      return out;
    }
    else {
      char* endptr = NULL;
      // Keep the pointer, because std::string doesn't necessarily
      // guarantee that this is the same across calls to c_str(), does
      // it?  Or perhaps it does...
      const char* t_ptr = t.c_str ();
      // We preset errno to 0, to distinguish success or failure after
      // calling strtoull.  Most implementations of the C standard
      // library written with threads in mind have errno be a macro that
      // expands to thread-local storage.  Thanks to the Linux
      // documentation for strtod ("man 3 strtod", Red Hat Enterprise
      // Linux 5) for advice with the following checks.
      errno = 0;
      const RealType val = rawConvert (t_ptr, &endptr);

      TEUCHOS_TEST_FOR_EXCEPTION(
        errno == ERANGE && (val != 0),
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<" << realTypeName
        << ", std::string>::convert: "
        "The value in the given string \"" << t << "\" overflows "
        << realTypeName << ".");
      //
      // mfh 20 Nov 2012: Should we treat underflow as an error?
      //
      TEUCHOS_TEST_FOR_EXCEPTION(
        errno == ERANGE && val == 0,
        std::invalid_argument,
        "Teuchos::ValueTypeConversionTraits<" << realTypeName
        << ", std::string>::convert: "
        "The value in the given string \"" << t << "\" underflows "
        << realTypeName << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        endptr == t_ptr, // See above discussion of c_str().
        std::invalid_argument,
        "Teuchos::ValueTypeConversionTraits<" << realTypeName
        << ", std::string>::convert: "
        "The conversion function was unable to read any floating-point data "
        "from the given string \"" << t << "\".");
      return val;
    }
  }

} // namespace (anonymous)


//
// Standard specializations of ValueTypeConversionTraits
//

//
// * Partial specialization for conversion from std::string to any type T.
//   There are full specializations for specific types T below.
//

/// \brief Convert an \c std::string to a type \c OutType.
///
/// This partial specialization assumes that \c OutType is default
/// constructible, and that <tt>operator>>(std::istream&,
/// OutType&)</tt> has been defined.  It does no bounds checking or
/// other input validation.
///
/// If you would like to add input validation for a specific output
/// type \c OutType, please implement a full specialization.  We
/// include many different full specializations in this file, so
/// please check this file first to see if we have already done the
/// work for you.
template<class OutType>
class ValueTypeConversionTraits<OutType, std::string> {
public:
  static OutType safeConvert (const std::string& t) {
    return convert (t);
  }

  static OutType convert (const std::string& t) {
    std::istringstream in (t);
    OutType out;
    t >> out;
    return out;
  }
};

//
// * Specializations for conversions from std::string to build-in
//   real-valued floating-point types.
//

//! Convert an \c std::string to a \c double.
template<>
class ValueTypeConversionTraits<double, std::string> {
public:
  static double convert (const std::string& t) {
    return realToString<double> (t, &strtod, "double");
  }

  static double safeConvert (const std::string& t) {
    return realToString<double> (t, &strtod, "double");
  }
};

//! Convert an \c std::string to a \c float.
template<>
class ValueTypeConversionTraits<float, std::string> {
public:
  static float convert (const std::string& t) {
#ifdef _ISOC99_SOURCE
    return realToString<float> (t, &strtof, "float");
#else
    // strtof is new in C99.  If you don't have it, just use strtod
    // and convert the resulting double to float.
    const double d = realToString<double> (t, &strtod, "double");
    return as<float> (d);
#endif // _ISOC99_SOURCE
  }

  static float safeConvert (const std::string& t) {
#ifdef _ISOC99_SOURCE
    return realToString<float> (t, &strtof, "float");
#else
    // strtof is new in C99.  If you don't have it, just use strtod
    // and convert the resulting double to float.
    const double d = realToString<double> (t, &strtod, "double");
    return asSafe<float> (d);
#endif // _ISOC99_SOURCE
  }
};

//! Convert an \c std::string to a <tt>long double</tt>.
template<>
class ValueTypeConversionTraits<long double, std::string> {
public:
  static long double convert (const std::string& t) {
#ifdef _ISOC99_SOURCE
    return realToString<long double> (t, &strtold, "long double");
#else
    // strtof is new in C99.  If you don't have it, just use
    // operator>>(std::istream&, long double&).
    return realToString<long double> (t, NULL, "long double");
#endif // _ISOC99_SOURCE
  }

  static long double safeConvert (const std::string& t) {
    return convert (t);
  }
};


//
// * Specializations for conversions from std::string to build-in integer types.
//

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

/// \brief Convert an \c std::string to a <tt>long long</tt>.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<long long, std::string> {
public:
  /// \brief Convert the given \c std::string to a <tt>long long</tt>, with checks.
  ///
  /// If the string overflows <tt>long long</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static long long safeConvert (const std::string& t) {
    return intToString<long long> (t, &strtoll, "long long");
  }

  //! Convert the given \c std::string to a <tt>long long</tt>.
  static long long convert (const std::string& t) {
    return safeConvert (t);
  }
};


/// \brief Convert an \c std::string to an <tt>unsigned long long</tt>.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<unsigned long long, std::string> {
public:
  /// \brief Convert the given \c std::string to an <tt>unsigned long long</tt>, with checks.
  ///
  /// If the string overflows <tt>unsigned long long</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static unsigned long long safeConvert (const std::string& t) {
    return intToString<unsigned long long> (t, &strtoull, "unsigned long long");
  }

  //! Convert the given \c std::string to an <tt>unsigned long long</tt>.
  static unsigned long long convert (const std::string& t) {
    return safeConvert (t);
  }
};

#endif // HAVE_TEUCHOS_LONG_LONG_INT


/// \brief Convert an \c std::string to a \c long.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<long, std::string> {
public:
  /// \brief Convert the given \c std::string to a \c long, with checks.
  ///
  /// If the string overflows <tt>long</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static long safeConvert (const std::string& t) {
    return intToString<long> (t, &strtol, "long");
  }

  //! Convert the given \c std::string to a \c long.
  static long convert (const std::string& t) {
    return safeConvert (t);
  }
};


/// \brief Convert an \c std::string to an <tt>unsigned long</tt>.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<unsigned long, std::string> {
public:
  /// \brief Convert the given std::string to an <tt>unsigned long</tt>, with checks.
  ///
  /// If the string overflows <tt>unsigned long</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static unsigned long safeConvert (const std::string& t) {
    return intToString<unsigned long> (t, &strtoul, "unsigned long");
  }

  //! Convert the given \c std::string to an <tt>unsigned long</tt>.
  static unsigned long convert (const std::string& t) {
    return safeConvert (t);
  }
};


/// \brief Convert an \c std::string to an \c int.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<int, std::string> {
private:
  //! Convert the given \c std::string to an intermediate \c long, with checks.
  static long safeConvertToLong (const std::string& t) {
    long val = 0;
    try {
      val = ValueTypeConversionTraits<long, std::string>::safeConvert (t);
    } catch (std::range_error&) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<int, std::string>::convert: "
        "The given std::string \"" << t << "\" is too big to fit into long, so there is no way it could fit into int.");
    } catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::invalid_argument,
        "Teuchos::ValueTypeConversionTraits<int, std::string>::convert: "
        "Intermediate conversion from std::string to long failed, with the following error message: "
        << e.what ());
    }
    return val;
  }

public:
  /// \brief Convert the given \c std::string to an <tt>int</tt>, with checks.
  ///
  /// If the string overflows <tt>int</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static int safeConvert (const std::string& t) {
    return asSafe<int> (safeConvertToLong (t));
  }

  //! Convert the given \c std::string to an \c int.
  static int convert (const std::string& t) {
    return as<int> (safeConvertToLong (t));
  }
};


/// \brief Convert an \c std::string to an <tt>unsigned int</tt>.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<unsigned int, std::string> {
private:
  //! Convert the given \c std::string to an intermediate <tt>unsigned long</tt>, with checks.
  static unsigned long safeConvertToUnsignedLong (const std::string& t) {
    unsigned long val = 0;
    try {
      val = as<unsigned long> (t);
    } catch (std::range_error&) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned int, std::string>::convert: "
        "The given std::string \"" << t << "\" is too big to fit into unsigned long, so there is no way it could fit into unsigned int.");
    } catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::invalid_argument,
        "Teuchos::ValueTypeConversionTraits<unsigned int, std::string>::convert: "
        "Intermediate conversion from std::string to unsigned long failed, with the following error message: "
        << e.what ());
    }
    return val;
  }

public:
  /// \brief Convert the given \c std::string to an <tt>unsigned int</tt>, with checks.
  ///
  /// If the string overflows <tt>unsigned int</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static unsigned int safeConvert (const std::string& t) {
    return asSafe<unsigned int> (safeConvertToUnsignedLong (t));
  }

  //! Convert the given \c std::string to an <tt>unsigned int</tt>.
  static unsigned int convert (const std::string& t) {
    return as<unsigned int> (safeConvertToUnsignedLong (t));
  }
};


/// \brief Convert an \c std::string to a \c short.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<short, std::string> {
private:
  //! Convert the given \c std::string to an intermediate \c long, with checks.
  static long safeConvertToLong (const std::string& t) {
    long val = 0;
    try {
      val = ValueTypeConversionTraits<long, std::string>::safeConvert (t);
    } catch (std::range_error&) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<short, std::string>::convert: "
        "The given std::string \"" << t << "\" is too big to fit into long, so there is no way it could fit into short.");
    } catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::invalid_argument,
        "Teuchos::ValueTypeConversionTraits<short, std::string>::convert: "
        "Intermediate conversion from std::string to long failed, with the following error message: "
        << e.what ());
    }
    return val;
  }

public:
  /// \brief Convert the given \c std::string to a <tt>short</tt>, with checks.
  ///
  /// If the string overflows <tt>short</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static short safeConvert (const std::string& t) {
    return asSafe<short> (safeConvertToLong (t));
  }

  //! Convert the given \c std::string to a \c short.
  static short convert (const std::string& t) {
    return as<int> (safeConvertToLong (t));
  }
};


/// \brief Convert an \c std::string to an <tt>unsigned short</tt>.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<unsigned short, std::string> {
private:
  //! Convert the given \c std::string to an intermediate <tt>unsigned long</tt>, with checks.
  static unsigned long safeConvertToUnsignedLong (const std::string& t) {
    unsigned long val = 0;
    try {
      val = as<unsigned long> (t);
    } catch (std::range_error&) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned short, std::string>::convert: "
        "The given std::string \"" << t << "\" is too big to fit into unsigned long, so there is no way it could fit into unsigned short.");
    } catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::invalid_argument,
        "Teuchos::ValueTypeConversionTraits<unsigned short, std::string>::convert: "
        "Intermediate conversion from std::string to unsigned long failed, with the following error message: "
        << e.what ());
    }
    return val;
  }

public:
  /// \brief Convert the given \c std::string to an <tt>unsigned short</tt>, with checks.
  ///
  /// If the string overflows <tt>unsigned short</tt>, this throws
  /// <tt>std::range_error</tt>.  If it does not contain an integer,
  /// this throws <tt>std::invalid_argument</tt>.
  static unsigned short safeConvert (const std::string& t) {
    return asSafe<unsigned short> (safeConvertToUnsignedLong (t));
  }

  //! Convert the given \c std::string to an <tt>unsigned short</tt>.
  static unsigned short convert (const std::string& t) {
    return as<unsigned short> (safeConvertToUnsignedLong (t));
  }
};

//
// * Specializations for conversions between built-in real-valued
//   floating-point types (like float and double).
//

//! Convert from \c double to \c float.
template<>
class ValueTypeConversionTraits<float, double> {
public:
  static float safeConvert (const double t) {
    // For floating-point types T, std::numeric_limits<T>::min()
    // returns the smallest positive value.  IEEE 754 types have a
    // sign bit, so the largest-magnitude negative value is the
    // negative of the largest-magnitude positive value.
    const float minVal = -std::numeric_limits<float>::max ();
    const float maxVal = std::numeric_limits<float>::max ();

    // NaN is neither less than nor greater than anything.  We just
    // let it pass through, per the rules for propagation of silent
    // NaN.  (Signaling NaN will signal, but that's OK.)
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<float, double>::safeConvert: "
      "Input double t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to float.");

    return static_cast<float> (t);
  }

  static float convert (const double t) {
    return static_cast<float> (t);
  }
};


//! Convert from <tt>long double</tt> to \c float.
template<>
class ValueTypeConversionTraits<float, long double> {
public:
  static float safeConvert (const long double t) {
    // For floating-point types T, std::numeric_limits<T>::min()
    // returns the smallest positive value.  IEEE 754 types have a
    // sign bit, so the largest-magnitude negative value is the
    // negative of the largest-magnitude positive value.
    const float minVal = -std::numeric_limits<float>::max ();
    const float maxVal = std::numeric_limits<float>::max ();

    // NaN is neither less than nor greater than anything.  We just
    // let it pass through, per the rules for propagation of silent
    // NaN.  (Signaling NaN will signal, but that's OK.)
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<float, long double>::safeConvert: "
      "Input long double t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to float.");

    return static_cast<float> (t);
  }

  static float convert (const long double t) {
    return static_cast<float> (t);
  }
};


//! Convert from <tt>long double</tt> to \c double.
template<>
class ValueTypeConversionTraits<double, long double> {
public:
  static double safeConvert (const long double t) {
    // For floating-point types T, std::numeric_limits<T>::min()
    // returns the smallest positive value.  IEEE 754 types have a
    // sign bit, so the largest-magnitude negative value is the
    // negative of the largest-magnitude positive value.
    const double minVal = -std::numeric_limits<double>::max ();
    const double maxVal = std::numeric_limits<double>::max ();

    // NaN is neither less than nor greater than anything.  We just
    // let it pass through, per the rules for propagation of silent
    // NaN.  (Signaling NaN will signal, but that's OK.)
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<double, long double>::safeConvert: "
      "Input long double t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to double.");

    return static_cast<double> (t);
  }

  static double convert (const long double t) {
    return static_cast<double> (t);
  }
};


//
// * Specializations for conversions from built-in real-valued
//   floating-point types (float and double) to build-in integer
//   types.
//

//! Convert from \c double to \c short.
template<>
class ValueTypeConversionTraits<short, double> {
public:
  /// \brief Convert the given \c double to a \c short.
  ///
  /// \warning Double-precision floating-point values may overflow
  ///   <tt>short</tt>.  You should use safeConvert() if you aren't sure
  ///   that the given value fits in an <tt>short</tt>.
  static short convert (const double t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<short> (t);
  }

  //! Convert the given \c double to a \c short, checking for overflow first.
  static short safeConvert (const double t) {
    const short minVal = std::numeric_limits<short>::min ();
    const short maxVal = std::numeric_limits<short>::max ();

    // Cases:
    // 1. sizeof(short) < sizeof(double) == 8
    // 2. sizeof(short) == sizeof(double) == 8
    // 3. sizeof(short) > sizeof(double) == 8
    //
    // Overflow when converting from double to short is possible only
    // for Case 1.  Loss of accuracy (rounding) is possible for Cases
    // 2 and 3, but safeConvert() only cares about overflow, not
    // rounding.  In Case 3, casting minVal or maxVal to double in
    // this case could result in overflow.  Thus, we only do the test
    // for Case 1.
    //
    // All three cases are legal according to both C++03 and C99.
    // However, I (mfh 15 Nov 2012) have never encountered Cases 2 and
    // 3.
    if (sizeof (short) < sizeof (double)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<short, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to short.");
    }
    return static_cast<short> (t);
  }
};


//! Convert from \c double to <tt>unsigned short</tt>.
template<>
class ValueTypeConversionTraits<unsigned short, double> {
public:
  //! Convert the given \c double to an <tt>unsigned short</tt>.
  static unsigned short convert (const double t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<unsigned short> (t);
  }

  //! Convert the given \c double to an <tt>unsigned short</tt>, checking for overflow first.
  static unsigned short safeConvert (const double t) {
    const unsigned short minVal = 0;
    const unsigned short maxVal = std::numeric_limits<unsigned short>::max ();

    if (sizeof (unsigned short) < sizeof (double)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned short, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned short.");
    }
    else { // Overflow isn't possible, but t < 0 isn't allowed.
      TEUCHOS_TEST_FOR_EXCEPTION(
       t < minVal,
       std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned short, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned short.");
    }
    return static_cast<unsigned short> (t);
  }
};


//! Convert from \c double to \c int.
template<>
class ValueTypeConversionTraits<int, double> {
public:
  /// \brief Convert the given \c double to an \c int.
  ///
  /// \warning Double-precision floating-point values may overflow
  ///   <tt>int</tt>.  You should use safeConvert() if you aren't sure
  ///   that the given value fits in an <tt>int</tt>.
  static int convert (const double t) {
    // Implicit conversion from double to int causes compiler
    // warnings, but static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert the given \c double to an \c int, checking for overflow first.
  static int safeConvert (const double t) {
    const int minVal = std::numeric_limits<int>::min ();
    const int maxVal = std::numeric_limits<int>::max ();

    // Cases:
    // 1. sizeof(int) < sizeof(double) == 8
    // 2. sizeof(int) == sizeof(double) == 8
    // 3. sizeof(int) > sizeof(double) == 8
    //
    // Overflow when converting from double to int is possible only
    // for Case 1.  Loss of accuracy (rounding) is possible for Cases
    // 2 and 3, but safeConvert() only cares about overflow, not
    // rounding.  Case 3 is quite rare, but casting minVal or maxVal
    // to double in this case could result in overflow.  Thus, we only
    // do the cast for Case 1.
    if (sizeof (int) < sizeof (double)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<int, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to int.");
    }
    return static_cast<int> (t);
  }
};


//! Convert from \c double to <tt>unsigned int</tt>.
template<>
class ValueTypeConversionTraits<unsigned int, double> {
public:
  //! Convert the given \c double to an <tt>unsigned int</tt>.
  static unsigned int convert (const double t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert the given \c double to an <tt>unsigned int</tt>, checking for overflow first.
  static unsigned int safeConvert (const double t) {
    const unsigned int minVal = 0;
    const unsigned int maxVal = std::numeric_limits<unsigned int>::max ();

    if (sizeof (unsigned int) < sizeof (double)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned int, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned int.");
    }
    else { // Overflow isn't possible, but t < 0 isn't allowed.
      TEUCHOS_TEST_FOR_EXCEPTION(
       t < minVal,
       std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned int, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned int.");
    }
    return static_cast<unsigned int> (t);
  }
};


//! Convert from \c double to \c long.
template<>
class ValueTypeConversionTraits<long, double> {
public:
  //! Convert the given \c double to \c long.
  static long convert (const double t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<long> (t);
  }

  //! Convert the given \c double to \c long, checking for overflow first.
  static long safeConvert (const double t) {
    const long minVal = std::numeric_limits<long>::min ();
    const long maxVal = std::numeric_limits<long>::max ();

    // Cases:
    // 1. sizeof(long) < sizeof(double) == 8
    // 2. sizeof(long) == sizeof(double) == 8
    // 3. sizeof(long) > sizeof(double) == 8
    //
    // Overflow when converting from double to long is possible only
    // for Case 1.  Loss of accuracy (rounding) is possible for Cases
    // 2 and 3, but safeConvert() only cares about overflow, not
    // rounding.  In Case 3, casting minVal or maxVal to double could
    // result in overflow.  Thus, we only test for Case 1.
    //
    // Case 1 is entirely possible, for example on Win64 (an
    // implementation of the LLP64 integer model, on which
    // sizeof(long) == 4, and sizeof(long long) == sizeof(void*) ==
    // 8).
    if (sizeof (long) < sizeof (double)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<long, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to long.");
    }
    return static_cast<long> (t);
  }
};


//! Convert from \c double to <tt>unsigned long</tt>.
template<>
class ValueTypeConversionTraits<unsigned long, double> {
public:
  //! Convert the given \c double to an <tt>unsigned long</tt>.
  static unsigned long convert (const double t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<unsigned long> (t);
  }

  //! Convert the given \c double to an <tt>unsigned long</tt>, checking for overflow first.
  static unsigned long safeConvert (const double t) {
    const unsigned long minVal = 0;
    const unsigned long maxVal = std::numeric_limits<unsigned long>::max ();

    if (sizeof (unsigned long) < sizeof (double)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned long, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned long.");
    }
    else { // Overflow isn't possible, but t < 0 isn't allowed.
      TEUCHOS_TEST_FOR_EXCEPTION(
       t < minVal,
       std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned long, double>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned long.");
    }
    return static_cast<unsigned long> (t);
  }
};

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

//! Convert from \c double to <tt>long long</tt>.
template<>
class ValueTypeConversionTraits<long long, double> {
public:
  //! Convert the given \c double to <tt>long long</tt>.
  static long long convert (const double t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<long long> (t);
  }

  //! Convert the given \c double to <tt>long long</tt>, checking for overflow first.
  static long safeConvert (const double t) {
    // Cases:
    // 1. sizeof(long long) < sizeof(double) == 8
    // 2. sizeof(long long) == sizeof(double) == 8
    // 3. sizeof(long long) > sizeof(double) == 8
    //
    // C99 (which defines long long) prohibits Case 1.  Case 2 could
    // result in loss of accuracy (rounding), but safeConvert() only
    // cares about overflow, not rounding.  In Case 3, casting minVal
    // or maxVal to double could result in overflow.  Thus, we don't
    // need to check anything here.
    return static_cast<long long> (t);
  }
};


//! Convert from \c double to <tt>unsigned long long</tt>.
template<>
class ValueTypeConversionTraits<unsigned long long, double> {
public:
  //! Convert the given \c double to <tt>unsigned long long</tt>.
  static unsigned long long convert (const double t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<unsigned long long> (t);
  }

  //! Convert the given \c double to <tt>unsigned long long</tt>, checking for overflow first.
  static long safeConvert (const double t) {
    // Cases:
    // 1. sizeof(unsigned long long) < sizeof(double) == 8
    // 2. sizeof(unsigned long long) == sizeof(double) == 8
    // 3. sizeof(unsigned long long) > sizeof(double) == 8
    //
    // C99 (which defines unsigned long long) prohibits Case 1.  Case
    // 2 could result in loss of accuracy (rounding), but
    // safeConvert() only cares about overflow, not rounding.  In Case
    // 3, casting minVal or maxVal to double could result in overflow.
    // Thus, we don't need to check the upper bound, though we still
    // need to check if the input is negative.
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < 0,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<unsigned long long, double>::safeConvert: "
      "Input double t = " << t << " is negative, which is invalid for conversion to unsigned long long.");
    return static_cast<unsigned long long> (t);
  }
};

#endif // HAVE_TEUCHOS_LONG_LONG_INT


//! Convert from \c float to \c short.
template<>
class ValueTypeConversionTraits<short, float> {
public:
  /// \brief Convert the given \c float to a \c short.
  ///
  /// \warning Single-precision floating-point values may overflow
  ///   <tt>short</tt>.  You should use safeConvert() if you aren't
  ///   sure that the given value fits in an <tt>short</tt>.
  static short convert (const float t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<short> (t);
  }

  //! Convert the given \c float to a \c short, checking for overflow first.
  static short safeConvert (const float t) {
    const short minVal = std::numeric_limits<short>::min ();
    const short maxVal = std::numeric_limits<short>::max ();

    // Cases:
    // 1. sizeof(short) < sizeof(float) == 4
    // 2. sizeof(short) == sizeof(float) == 4
    // 3. sizeof(short) > sizeof(float) == 4
    //
    // Overflow when converting from float to short is possible only
    // for Case 1.  Loss of accuracy (rounding) is possible for Cases
    // 2 and 3, but safeConvert() only cares about overflow, not
    // rounding.  In Case 3, casting minVal or maxVal to float in this
    // case could result in overflow.  Thus, we only do the test for
    // Case 1.
    //
    // All three cases are legal according to both C++03 and C99.  I
    // (mfh 15 Nov 2012) think Case 1 is the most common, but Case 2
    // is certainly reasonable.  (For example, some hardware prefers
    // to work only with 32-bit words, so _every_ built-in type has
    // size a multiple of 4 bytes.)
    if (sizeof (short) < sizeof (float)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<short, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to short.");
    }

    return static_cast<short> (t);
  }
};


//! Convert from \c float to <tt>unsigned short</tt>.
template<>
class ValueTypeConversionTraits<unsigned short, float> {
public:
  //! Convert the given \c float to an <tt>unsigned short</tt>.
  static unsigned short convert (const float t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<unsigned short> (t);
  }

  //! Convert the given \c float to an <tt>unsigned short</tt>, checking for overflow first.
  static unsigned short safeConvert (const float t) {
    const unsigned short minVal = 0;
    const unsigned short maxVal = std::numeric_limits<unsigned short>::max ();

    if (sizeof (unsigned short) < sizeof (float)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned short, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned short.");
    }
    else { // Overflow isn't possible, but t < 0 isn't allowed.
      TEUCHOS_TEST_FOR_EXCEPTION(
       t < static_cast<float> (minVal),
       std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned short, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned short.");
    }
    return static_cast<unsigned short> (t);
  }
};


//! Convert from \c float to \c int.
template<>
class ValueTypeConversionTraits<int, float> {
public:
  //! Convert the given \c float to an \c int.
  static int convert (const float t) {
    // Implicit conversion from float to int may cause compiler
    // warnings, but static_cast does not.  Overflow here would mean
    // that sizeof(int) < sizeof(float), which is legal, but unlikely
    // on platforms of interest.
    return static_cast<int> (t);
  }

  //! Convert the given \c float to an \c int.
  static int safeConvert (const float t) {
    const int minVal = std::numeric_limits<int>::min ();
    const int maxVal = std::numeric_limits<int>::max ();

    // Cases:
    // 1. sizeof(int) < sizeof(float) == 4
    // 2. sizeof(int) == sizeof(float) == 4
    // 3. sizeof(int) > sizeof(float) == 4
    //
    // Overflow when converting from float to int is possible only for
    // Case 1.  Loss of accuracy (rounding) is possible for Cases 2
    // and 3, but safeConvert() only cares about overflow, not
    // rounding.  Case 3 is rare, but casting minVal or maxVal to
    // float in this case could result in loss of accuracy
    // (sizeof(int) == 8 or 16) or overflow (sizeof(int) > 16).  Thus,
    // we only do the test for Case 1.
    if (sizeof (unsigned int) < sizeof (float)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<int, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range ["
        << minVal << ", " << maxVal << "] for conversion to int.");
    }
    return static_cast<int> (t);
  }
};


//! Convert from \c float to <tt>unsigned int</tt>.
template<>
class ValueTypeConversionTraits<unsigned int, float> {
public:
  //! Convert the given \c float to an <tt>unsigned int</tt>.
  static unsigned int convert (const float t) {
    // Implicit conversion from float to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert the given \c float to an <tt>unsigned int</tt>, checking first or under- or overflow.
  static unsigned int safeConvert (const float t) {
    const unsigned int minVal = 0; // Had better be, since it's unsigned.
    const unsigned int maxVal = std::numeric_limits<unsigned int>::max ();

    if (sizeof (unsigned int) < sizeof (float)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned int, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned int.");
    }
    else { // Overflow isn't possible, but t < 0 isn't allowed.
      TEUCHOS_TEST_FOR_EXCEPTION(
       t < minVal,
       std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned int, float>::safeConvert: "
        "Input double t = " << t << " is out of the valid range [" << minVal
        << ", " << maxVal << "] for conversion to unsigned int.");
    }
    return static_cast<unsigned int> (t);
  }
};


//! Convert from \c float to \c long.
template<>
class ValueTypeConversionTraits<long, float> {
public:
  //! Convert the given \c float to an \c long.
  static long convert (const float t) {
    // Implicit conversion from float to long may cause compiler
    // warnings, but static_cast does not.  Overflow here would mean
    // that sizeof(long) < sizeof(float), which is legal, but unlikely
    // on platforms of longerest.
    return static_cast<long> (t);
  }

  //! Convert the given \c float to an \c long, checking first for overflow.
  static long safeConvert (const float t) {
    const long minVal = std::numeric_limits<long>::min ();
    const long maxVal = std::numeric_limits<long>::max ();

    // Cases:
    // 1. sizeof(long) < sizeof(float) == 4
    // 2. sizeof(long) == sizeof(float) == 4
    // 3. sizeof(long) > sizeof(float) == 4
    //
    // Overflow when converting from float to long is possible only
    // for Case 1.  Loss of accuracy (rounding) is possible for Cases
    // 2 and 3, but safeConvert() only cares about overflow, not
    // rounding.  Casting minVal or maxVal to double in Case 3 could
    // result in overflow.  Thus, we only do the cast for Case 1.
    //
    // I've never encountered a Case 1 platform (mfh 14 Nov 2012).
    // C99 actually forbids it, though I don't think a valid C++
    // compiler (for version C++03 of the language standard) needs to
    // implement C99 (mfh 14 Nov 2012).  Case 2 occurs in Win64
    // (64-bit Windows) and other implementations of (I32L32)LLP64.
    // Case 3 is common (e.g., in the (I32)LP64 integer model of
    // GNU/Linux and other operating systems).
    if (sizeof (long) < sizeof (float)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < minVal || t > maxVal,
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<long, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range ["
        << minVal << ", " << maxVal << "] for conversion to long.");
    }
    return static_cast<long> (t);
  }
};


//! Convert from \c float to <tt>unsigned long</tt>.
template<>
class ValueTypeConversionTraits<unsigned long, float> {
public:
  //! Convert the given \c float to an <tt>unsigned long</tt>.
  static unsigned long convert (const float t) {
    // Implicit conversion from float to unsigned long may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned long> (t);
  }

  //! Convert the given \c float to an <tt>unsigned long</tt>, checking first or under- or overflow.
  static unsigned long safeConvert (const float t) {
    const unsigned long minVal = 0; // Had better be, since it's unsigned.
    const unsigned long maxVal = std::numeric_limits<unsigned long>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<unsigned long, float>::safeConvert: "
      << std::endl
      << "Input float t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to unsigned long.  "
      << std:: endl
      << "(Remember that unsigned types cannot represent negative values.)");

    return static_cast<unsigned long> (t);
  }
};

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

//! Convert from \c float to <tt>long long</tt>.
template<>
class ValueTypeConversionTraits<long long, float> {
public:
  //! Convert the given \c float to a <tt>long long</tt>.
  static long long convert (const float t) {
    return static_cast<long long> (t);
  }

  //! Convert the given \c float to a <tt>long long</tt>, checking first for overflow.
  static long long safeConvert (const float t) {
    // The C99 standard (Section 5.2.4.2.1) actually requires
    // sizeof(long long) >= 64, so overflow is impossible.
    return static_cast<long long> (t);
  }
};


//! Convert from \c float to <tt>unsigned long long</tt>.
template<>
class ValueTypeConversionTraits<unsigned long long, float> {
public:
  //! Convert the given \c float to an <tt>unsigned long long</tt>.
  static unsigned long long convert (const float t) {
    return static_cast<unsigned long long> (t);
  }

  //! Convert the given \c float to an <tt>unsigned long long</tt>, checking first for overflow.
  static unsigned long long safeConvert (const float t) {
    // The C99 standard (Section 5.2.4.2.1) actually requires
    // sizeof(long long) >= 64, so overflow is impossible.  However,
    // we still forbid negative inputs, since the target type is
    // unsigned.
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < 0,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<unsigned long long, float>::safeConvert: "
      "Input float t = " << t << " is negative, which is invalid for conversion to unsigned long long.");
    return static_cast<unsigned long long> (t);
  }
};

#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// * Specializations for conversions between a unsigned built-in
//   integer type and the signed version of the same type (either
//   direction).
//

namespace {
// Implementation of conversion from an unsigned built-in integer
// type, to an signed built-in integer type with the same number of
// bits.
template<class SignedIntType, class UnsignedIntType>
class UnsignedToSignedValueTypeConversionTraits {
public:
  /// Convert from unsigned to signed.
  ///
  /// \warning Some unsigned integer values may overflow the signed
  ///   integer type, resulting in a negative number when the original
  ///   number was positive.  You should use safeConvert() if you
  ///   aren't sure that the given unsigned value fits in the signed
  ///   type.
  static SignedIntType convert (const UnsignedIntType t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<SignedIntType> (t);
  }

  //! Convert from unsigned to signed, checking for overflow first.
  static SignedIntType safeConvert (const UnsignedIntType t) {
    using Teuchos::TypeNameTraits;
    const SignedIntType maxSigned = std::numeric_limits<SignedIntType>::max ();

    // SignedIntType and UnsignedIntType have the same number of bits,
    // so it suffices (via two's complement arithmetic) to check
    // whether the cast turned a positive number negative.
    const SignedIntType signedVal = static_cast<SignedIntType> (t);
    TEUCHOS_TEST_FOR_EXCEPTION(
      signedVal < 0,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<" << TypeNameTraits<SignedIntType>::name ()
      << ", " << TypeNameTraits<UnsignedIntType>::name () << ">::safeConvert: "
      "Input " << TypeNameTraits<UnsignedIntType>::name () << " t = " << t
      << " is out of the valid range [0, " << ", " << maxSigned
      << "] for conversion to " << TypeNameTraits<SignedIntType>::name () << ".");
    return signedVal;
  }
};


// Implementation of conversion from a signed built-in integer type,
// to an unsigned built-in integer type with the same number of bits.
template<class UnsignedIntType, class SignedIntType>
class SignedToUnsignedValueTypeConversionTraits {
public:
  //! Convert the given unsigned integer to a signed integer of the same size.
  static UnsignedIntType convert (const SignedIntType t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<UnsignedIntType> (t);
  }

  //! Convert from signed to unsigned, checking for underflow first.
  static UnsignedIntType safeConvert (const SignedIntType t) {
    using Teuchos::TypeNameTraits;

    // SignedIntType and UnsignedIntType have the same number of bits,
    // so it suffices (via two's complement arithmetic) to check
    // whether the input is negative.
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<SignedIntType> (0),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<" << TypeNameTraits<UnsignedIntType>::name ()
      << ", " << TypeNameTraits<SignedIntType>::name () << ">::safeConvert: "
      "Input " << TypeNameTraits<SignedIntType>::name () << " t = " << t
      << " is negative, so it cannot be correctly converted to the unsigned type "
      << TypeNameTraits<UnsignedIntType>::name () << ".");

    return static_cast<UnsignedIntType> (t);
  }
};

} // namespace (anonymous)


//! Convert from <tt>unsigned short<tt> to \c short.
template<>
class ValueTypeConversionTraits<short, unsigned short> {
public:
  static short convert (const unsigned short t) {
    return UnsignedToSignedValueTypeConversionTraits<short, unsigned short>::convert (t);
  }

  static short safeConvert (const unsigned short t) {
    return UnsignedToSignedValueTypeConversionTraits<short, unsigned short>::safeConvert (t);
  }
};


//! Convert from <tt>short<tt> to <tt>unsigned short</tt>.
template<>
class ValueTypeConversionTraits<unsigned short, short> {
public:
  static unsigned short convert (const short t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned short, short>::convert (t);
  }

  static unsigned short safeConvert (const short t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned short, short>::safeConvert (t);
  }
};


//! Convert from <tt>unsigned int<tt> to \c int.
template<>
class ValueTypeConversionTraits<int, unsigned int> {
public:
  static int convert (const unsigned int t) {
    return UnsignedToSignedValueTypeConversionTraits<int, unsigned int>::convert (t);
  }

  static int safeConvert (const unsigned int t) {
    return UnsignedToSignedValueTypeConversionTraits<int, unsigned int>::safeConvert (t);
  }
};


//! Convert from <tt>int<tt> to <tt>unsigned int</tt>.
template<>
class ValueTypeConversionTraits<unsigned int, int> {
public:
  static unsigned int convert (const int t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned int, int>::convert (t);
  }

  static unsigned int safeConvert (const int t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned int, int>::safeConvert (t);
  }
};


//! Convert from <tt>unsigned long<tt> to \c long.
template<>
class ValueTypeConversionTraits<long, unsigned long> {
public:
  static long convert (const unsigned long t) {
    return UnsignedToSignedValueTypeConversionTraits<long, unsigned long>::convert (t);
  }

  static long safeConvert (const unsigned long t) {
    return UnsignedToSignedValueTypeConversionTraits<long, unsigned long>::safeConvert (t);
  }
};


//! Convert from <tt>long<tt> to <tt>unsigned long</tt>.
template<>
class ValueTypeConversionTraits<unsigned long, long> {
public:
  static unsigned long convert (const long t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned long, long>::convert (t);
  }

  static unsigned long safeConvert (const long t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned long, long>::safeConvert (t);
  }
};


#ifdef HAVE_TEUCHOS_LONG_LONG_INT

//! Convert from <tt>unsigned long long<tt> to <tt>long long</tt>.
template<>
class ValueTypeConversionTraits<long long, unsigned long long> {
public:
  static long long convert (const unsigned long long t) {
    return UnsignedToSignedValueTypeConversionTraits<long long, unsigned long long>::convert (t);
  }

  static long long safeConvert (const unsigned long long t) {
    return UnsignedToSignedValueTypeConversionTraits<long long, unsigned long long>::safeConvert (t);
  }
};


//! Convert from <tt>long long<tt> to <tt>unsigned long long</tt>.
template<>
class ValueTypeConversionTraits<unsigned long long, long long> {
public:
  static unsigned long long convert (const long long t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned long long, long long>::convert (t);
  }

  static unsigned long long safeConvert (const long long t) {
    return SignedToUnsignedValueTypeConversionTraits<unsigned long long, long long>::safeConvert (t);
  }
};

#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// * Specializations for conversions between different built-in
//   integer types.
//

//! Convert from \c int to \c short.
template<>
class ValueTypeConversionTraits<short, int> {
public:
  /// \brief Convert the given \c int to a \c short.
  ///
  /// \warning \c int values may overflow \c short, depending on your
  ///   platform.  You should use safeConvert() if you aren't sure
  ///   that the given \c int value fits in a \c short.
  static short convert (const int t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<short> (t);
  }

  //! Convert from \c int to \c short, checking for overflow first.
  static short safeConvert (const int t) {
    const short minShort = std::numeric_limits<short>::min ();
    const short maxShort = std::numeric_limits<short>::max ();

    // Casting from short to int never overflows, since the C++
    // standard guarantees that sizeof (short) <= sizeof (int).
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<int> (minShort) ||
      t > static_cast<int> (maxShort),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<short, int>::safeConvert: "
      "Input int t = " << t << " is out of the valid range [" << minShort
      << ", " << maxShort << "] for conversion to short.");

    return static_cast<short> (t);
  }
};


//! Convert from \c long to \c short.
template<>
class ValueTypeConversionTraits<short, long> {
public:
  /// \brief Convert the given \c long to a \c short.
  ///
  /// \warning \c long integer values may overflow \c short, depending
  ///   on your platform.  You should use safeConvert() if you aren't
  ///   sure that the given \c long value fits in a \c short.
  static short convert (const long t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<short> (t);
  }

  //! Convert from \c long to \c short, checking for overflow first.
  static short safeConvert (const long t) {
    const short minShort = std::numeric_limits<short>::min ();
    const short maxShort = std::numeric_limits<short>::max ();

    // Casting from short to long never overflows, since the C++
    // standard guarantees that sizeof (short) <= sizeof (long).
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<long> (minShort) ||
      t > static_cast<long> (maxShort),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<short, long>::safeConvert: "
      "Input long t = " << t << " is out of the valid range [" << minShort
      << ", " << maxShort << "] for conversion to short.");

    return static_cast<short> (t);
  }
};


//! Convert from \c long to \c int.
template<>
class ValueTypeConversionTraits<int, long> {
public:
  /// \brief Convert the given \c long to an \c int.
  ///
  /// \warning \c long integer values may overflow \c int, depending
  ///   on your platform.  You should use safeConvert() if you aren't
  ///   sure that the given \c long value fits in an \c int.
  static int convert (const long t) {
    // Implicit conversion from long to int may cause compiler
    // warnings, but static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert from \c long to \c int, checking for overflow first.
  static int safeConvert (const long t) {
    const int minInt = std::numeric_limits<int>::min ();
    const int maxInt = std::numeric_limits<int>::max ();

    // Casting from int to long never overflows, since the C++
    // standard guarantees that sizeof (int) <= sizeof (long).
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<long> (minInt) ||
      t > static_cast<long> (maxInt),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<int, long>::safeConvert: "
      "Input long t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to int.");

    // Implicit conversion from long to int may cause compiler
    // warnings, but static_cast does not.
    return static_cast<int> (t);
  }
};


//! Convert from <tt>unsigned long</tt> to \c int.
template<>
class ValueTypeConversionTraits<int, unsigned long> {
public:
  /// \brief Convert the given <tt>unsigned long</tt> to an \c int.
  ///
  /// \warning <tt>unsigned long</tt> values may overflow
  ///   <tt>int</tt>, depending on your platform.  You should use
  ///   safeConvert() if you aren't sure that the given <tt>unsigned
  ///   long</tt> value fits in an <tt>int</tt>.
  static int convert (const unsigned long t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert from <tt>unsigned long</tt> to \c int, checking for overflow first.
  static int safeConvert (const unsigned long t) {
    const int minInt = std::numeric_limits<int>::min ();
    const int maxInt = std::numeric_limits<int>::max ();

    // On some platforms, sizeof(int) == sizeof(long).  (This is the
    // "LLP64" model of Win64, which aims for backwards compatibility
    // with 32-bit code by making sizeof(int) == sizeof(long) == 4.)
    // If this is the case, then we can't safely cast unsigned long to
    // int, or unsigned int to long, because values with the most
    // significant bit set will overflow to negative values.

    // The C++ standard promises that sizeof (int) <= sizeof (unsigned
    // long).  We use #if with INT_MAX and LONG_MAX to test for this,
    // rather than if statements, in order to avoid a compiler
    // warning.  Thanks to Jeremie Gaidamour (13 Nov 2012) for letting
    // me know about the warning.
#if INT_MAX == LONG_MAX
    // The two types have the same number of bits.  Thus,
    // two's-complement arithmetic means that if casting from unsigned
    // long to int results in a negative number, it overflowed.
    // Otherwise, it didn't overflow (same number of bits).
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<int> (t) < static_cast<int> (0),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<int, unsigned long>::safeConvert: "
      "Input unsigned long t = " << t << " is out of the valid range ["
      << minInt << ", " << maxInt << "] for conversion to int.");
#else // INT_MAX < LONG_MAX
    // t is unsigned, so it is >= 0 by definition.
    // Casting from int to unsigned long won't overflow in this case.
    TEUCHOS_TEST_FOR_EXCEPTION(
      t > static_cast<unsigned long> (maxInt),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<int, unsigned long>::safeConvert: "
      "Input unsigned long t = " << t << " is out of the valid range ["
      << minInt << ", " << maxInt << "] for conversion to int.  An unchecked "
      "cast would have resulted in " << static_cast<int> (t) << ".");
#endif // INT_MAX == LONG_MAX

    // Implicit conversion from unsigned long to int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<int> (t);
  }
};


//! Convert from <tt>unsigned int</tt> to <tt>long</tt>.
template<>
class ValueTypeConversionTraits<long, unsigned int> {
public:
  /// \brief Convert the given <tt>unsigned int</tt> to a \c long.
  ///
  /// \warning On some platforms (e.g., Windows, or any other platform
  ///   that implements the LLP64 model), <tt>unsigned int</tt>
  ///   integer values may overflow <tt>long</tt>.  You should use
  ///   safeConvert() if you aren't sure that the given <tt>unsigned
  ///   int</tt> value fits in a <tt>long</tt>.
  static long convert (const unsigned int t) {
    // Implicit conversion may cause compiler warnings, but
    // static_cast does not.
    return static_cast<long> (t);
  }

  //! Convert from <tt>unsigned int</tt> to \c long, checking for overflow first.
  static long safeConvert (const unsigned int t) {
    // On some platforms, sizeof(int) == sizeof(long).  (This is the
    // "LLP64" model of Win64, which aims for backwards compatibility
    // with 32-bit code by making sizeof(int) == sizeof(long) == 4.)
    // If this is the case, then we can't safely cast unsigned long to
    // int, or unsigned int to long, because values with the most
    // significant bit set will overflow to negative values.

    // The C++ standard promises that sizeof (unsigned int) <= sizeof
    // (long).  If strictly less, then the conversion won't overflow.
    // We protect the test with an #ifdef ... #endif to avoid compiler
    // warnings like the following: "warning: comparison is always
    // false due to limited range of data type".
#if UINT_MIN == LONG_MIN && UINT_MAX == LONG_MAX
    const long minLong = std::numeric_limits<long>::min ();
    const long maxLong = std::numeric_limits<long>::max ();

    // The two types have the same number of bits.  Thus,
    // two's-complement arithmetic means that if casting from
    // unsigned int to long results in a negative number, it
    // overflowed.  Otherwise, it didn't overflow (same number of
    // bits).
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<long> (t) < static_cast<long> (0),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<long, unsigned int>::safeConvert: "
      "Input unsigned int t = " << t << " is out of the valid range ["
      << minLong << ", " << maxLong << "] for conversion to long.");
#endif // UINT_MAX == LONG_MAX

    return static_cast<long> (t);
  }
};


//! Convert from \c long to <tt>unsigned int</tt>.
template<>
class ValueTypeConversionTraits<unsigned int, long> {
public:
  /// \brief Convert the given \c long to an <tt>unsigned int</tt>.
  ///
  /// \warning \c long integer values may overflow <tt>unsigned
  ///   int</tt>, depending on your platform.  You should use
  ///   safeConvert() if you aren't sure that the given \c long value
  ///   fits in an <tt>unsigned int</tt>.
  static unsigned int convert (const long t) {
    // Implicit conversion from long to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert from \c long to <tt>unsigned int</tt>, checking for underflow or overflow first.
  static unsigned int safeConvert (const long t) {
    // On some platforms, sizeof(int) == sizeof(long).  (This is the
    // "LLP64" model of Win64, which aims for backwards compatibility
    // with 32-bit code by making sizeof(int) == sizeof(long) == 4.)
    // In this case, conversion from long to unsigned int can't
    // overflow.

    // The C++ standard promises that sizeof (unsigned int) <= sizeof (long).
    if (sizeof (unsigned int) < sizeof (long)) {
      const unsigned int maxInt = std::numeric_limits<unsigned int>::max ();

      TEUCHOS_TEST_FOR_EXCEPTION(
        t < static_cast<long> (0) || t > static_cast<long> (maxInt),
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned int, long>::safeConvert: "
        "Input long t = " << t << " is out of the valid range [0, "
        << maxInt << "] for conversion to unsigned int.");
    }
    // Implicit conversion from long to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }
};


//! Convert from <tt>unsigned long</tt> to <tt>unsigned int</tt>.
template<>
class ValueTypeConversionTraits<unsigned int, unsigned long> {
public:
  /// \brief Convert the given <tt>unsigned long</tt> to an <tt>unsigned int</tt>.
  ///
  /// \warning <tt>unsigned long</tt> integer values may overflow
  ///   <tt>unsigned int</tt>, depending on your platform.  You should
  ///   use safeConvert() if you aren't sure that the given
  ///   <tt>unsigned long</tt> value fits in an <tt>unsigned int</tt>.
  static unsigned int convert (const unsigned long t) {
    // Implicit conversion from unsigned long to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert from <tt>unsigned long</tt> to <tt>unsigned int</tt>, checking for overflow first.
  static unsigned int safeConvert (const unsigned long t) {
    const unsigned int minInt = 0; // Had better be, since it's unsigned.
    const unsigned int maxInt = std::numeric_limits<unsigned int>::max ();

    // t >= 0 by definition, because it is unsigned.
    TEUCHOS_TEST_FOR_EXCEPTION(
      t > static_cast<unsigned long> (maxInt),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<unsigned int, unsigned long>::safeConvert: "
      "Input unsigned long t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to unsigned int.");

    // Implicit conversion from unsigned long to unsigned int may
    // cause compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }
};

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

//! Convert from <tt>long long</tt> to \c int.
template<>
class ValueTypeConversionTraits<int, long long> {
public:
  /// \brief Convert the given <tt>long long</tt> to an \c int.
  ///
  /// \warning <tt>long long</tt> integer values may overflow \c int.
  ///   You should use safeConvert() if you aren't sure that the given
  ///   value fits in an \c int.
  static int convert (const long long t) {
    // Implicit conversion from long long to int may cause compiler
    // warnings, but static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert from <tt>long long</tt> to <tt>int</tt>, checking for overflow first.
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


//! Convert from <tt>long long</tt> to <tt>unsigned int</tt>.
template<>
class ValueTypeConversionTraits<unsigned int, long long> {
public:
  /// \brief Convert the given <tt>long long</tt> to an </tt>unsigned int</tt>.
  ///
  /// \warning <tt>long long</tt> integer values may overflow
  ///   </tt>unsigned int</tt>.  You should use safeConvert() if you
  ///   aren't sure that the given value fits in an <tt>unsigned
  ///   int</tt>.
  static unsigned int convert (const long long t) {
    // Implicit conversion from long long to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert from <tt>long long</tt> to <tt>unsigned int</tt>, checking for overflow first.
  static unsigned int safeConvert (const long long t) {
    const unsigned int minInt = 0; // Had better be, because it's unsigned.
    const unsigned int maxInt = std::numeric_limits<unsigned int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < static_cast<long long> (minInt) || t > static_cast<long long> (maxInt),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<unsigned int, long long>::safeConvert: "
      "Input long long t = " << t << " is out of the valid range [" << minInt
      << ", " << maxInt << "] for conversion to unsigned int.");

    // Implicit conversion from long long to unsigned int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }
};


//! Convert from <tt>unsigned long long</tt> to <tt>int</tt>.
template<>
class ValueTypeConversionTraits<int, unsigned long long> {
public:
  /// \brief Convert the given <tt>unsigned long long</tt> to an <tt>int</tt>.
  ///
  /// \warning <tt>unsigned long long</tt> integer values may overflow
  ///   \c int.  You should use safeConvert() if you aren't sure that
  ///   the given value fits in an \c int.
  static int convert (const unsigned long long t) {
    // Implicit conversion from unsigned long long to int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<int> (t);
  }

  //! Convert from <tt>unsigned long long</tt> to <tt>int</tt>, checking for overflow first.
  static int safeConvert (const unsigned long long t) {
    const int minInt = std::numeric_limits<int>::min ();
    const int maxInt = std::numeric_limits<int>::max ();

    // t >= 0 by definition, because it is unsigned.
    TEUCHOS_TEST_FOR_EXCEPTION(
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


//! Convert from <tt>unsigned long long</tt> to </tt>unsigned int</tt>.
template<>
class ValueTypeConversionTraits<unsigned int, unsigned long long> {
public:
  /// \brief Convert the given <tt>unsigned long long</tt> to an <tt>unsigned int</tt>.
  ///
  /// \warning <tt>unsigned long long</tt> integer values may overflow
  ///   <tt>unsigned int</tt>.  You should use safeConvert() if you
  ///   aren't sure that the given value fits in an <tt>unsigned
  ///   int</tt>.
  static unsigned int convert (const unsigned long long t) {
    // Implicit conversion from unsigned long long to unsigned int may
    // cause compiler warnings, but static_cast does not.
    return static_cast<unsigned int> (t);
  }

  //! Convert from <tt>unsigned long long</tt> to <tt>unsigned int</tt>, checking for overflow first.
  static unsigned int safeConvert (const unsigned long long t) {
    const unsigned int minInt = 0; // Had better be, since it's unsigned.
    const unsigned int maxInt = std::numeric_limits<unsigned int>::max ();

    // t >= 0 by definition, because it is unsigned.
    TEUCHOS_TEST_FOR_EXCEPTION(
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

//
// * Conversions from built-in integer types to built-in real-valued
//   floating-point types.
//

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

//! Convert from <tt>long long</tt> to \c float.
template<>
class ValueTypeConversionTraits<float, long long> {
public:
  /// \brief Convert the given <tt>long long</tt> to a \c float.
  ///
  /// \warning <tt>long long</tt> integer values may overflow
  ///   <tt>float</tt>.  You should use safeConvert() if you aren't
  ///   sure that the given value fits in a <tt>float</tt>.
  static float convert (const long long t) {
    // Implicit conversion from long long to float may cause compiler
    // warnings, but static_cast does not.
    return static_cast<float> (t);
  }

  //! Convert from <tt>long long</tt> to \c float, checking for overflow first.
  static float safeConvert (const long long t) {
    // std::numeric_limits<float>::min() gives the minimum _positive_
    // normalized value of type float.  IEEE 754 floating-point values
    // can change sign just by flipping the sign bit, so the "most
    // negative" finite float is just the negative of the "most
    // positive" finite float.
    const float minFloat = -std::numeric_limits<float>::max ();
    const float maxFloat = std::numeric_limits<float>::max ();

    // mfh 16 Nov 2012: On my platform (gcc 4.7.2, Red Hat Linux 5,
    // Intel x86_64), first casting [minFloat,maxFloat] to long long
    // (so that the comparison only compares long long values)
    // gives different results in the comparison below than just
    // comparing t (as a long long) with minFloat and maxFloat.  It
    // doesn't matter whether you use static_cast<long long> (...) or
    // (long long) (...) to do the cast: the original float interval
    // of [-3.40282e+38, 3.40282e+38] becomes [-9223372036854775808,
    // -9223372036854775808], which is obviously wrong.
    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minFloat || t > maxFloat,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<float, long long>::safeConvert: "
      "Input long long t = " << t << " is out of the valid range [" << minFloat
      << ", " << maxFloat << "] for conversion to float.");

    // Implicit conversion from long long to float may cause compiler
    // warnings, but static_cast does not.
    return static_cast<float> (t);
  }
};


//! Convert from <tt>unsigned long long</tt> to \c float.
template<>
class ValueTypeConversionTraits<float, unsigned long long> {
public:
  /// \brief Convert the given <tt>unsigned long long</tt> to a \c float.
  ///
  /// \warning <tt>unsigned long long</tt> integer values may overflow
  ///   \c float.  You should use safeConvert() if you aren't sure
  ///   that the given value fits in a \c float.
  static float convert (const unsigned long long t) {
    // Implicit conversion from unsigned long long to float may cause
    // compiler warnings, but static_cast does not.
    return static_cast<float> (t);
  }

  //! Convert from <tt>unsigned long long</tt> to \c float, checking for overflow first.
  static float safeConvert (const unsigned long long t) {
    // std::numeric_limits<float>::min() gives the minimum _positive_
    // normalized value of type float.  IEEE 754 floating-point values
    // can change sign just by flipping the sign bit, so the "most
    // negative" finite float is just the negative of the "most
    // positive" finite float.
    const float minFloat = -std::numeric_limits<float>::max ();
    const float maxFloat = std::numeric_limits<float>::max ();

    // t >= 0 by definition, because it is unsigned.
    //
    // mfh 16 Nov 2012: See my note above on the <float, long long>
    // specialization that explains why I don't cast maxFloat to
    // unsigned long long here.
    TEUCHOS_TEST_FOR_EXCEPTION(
      t > maxFloat,
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<float, unsigned long long>::safeConvert: "
      "Input unsigned long long t = " << t << " is out of the valid range [" << minFloat
      << ", " << maxFloat << "] for conversion to float.");

    // Implicit conversion from unsigned long long to float may cause
    // compiler warnings, but static_cast does not.
    return static_cast<float> (t);
  }
};

#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// * Other conversions
//

/** \brief Convert raw C string to std::string. */
template<int N>
class ValueTypeConversionTraits<std::string, char[N]> {
public:
  static std::string convert( const char t[] )
    { return std::string(t); }
  static std::string safeConvert( const char t[] )
    { return std::string(t); }
};

//
// * Conversions for dd_real and qd_real
//

#ifdef HAVE_TEUCHOS_QD

/** \brief Convert qd_real to double. */
template <>
class ValueTypeConversionTraits<double, qd_real> {
public:
  inline static double convert (const qd_real t) {
    return to_double (t);
  }
  static double safeConvert (const qd_real t) {
    // std::numeric_limits<double>::min() gives the minimum _positive_
    // normalized value of type double.  IEEE 754 floating-point
    // values can change sign just by flipping the sign bit, so the
    // "most negative" finite double is just the negative of the "most
    // positive" finite double.
    const qd_real minVal = -std::numeric_limits<double>::max ();
    const qd_real maxVal = std::numeric_limits<double>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<double, qd_real>::safeConvert: "
      "Input qd_real t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to double.");

    return to_double (t);
  }
};

/** \brief Convert qd_real to float. */
template <>
class ValueTypeConversionTraits<float, qd_real> {
public:
  inline static float convert (const qd_real t) {
    // In a debug build, this should also test the double->float
    // conversion for overflow.
    return as<float> (to_double (t));
  }

  static float safeConvert (const qd_real t) {
    // std::numeric_limits<float>::min() gives the minimum _positive_
    // normalized value of type float.  IEEE 754 floating-point
    // values can change sign just by flipping the sign bit, so the
    // "most negative" finite float is just the negative of the "most
    // positive" finite float.
    //
    // qd_real has a constructor for double, but not for float,
    // so we cast to double first.
    const qd_real minVal = static_cast<double> (-std::numeric_limits<float>::max ());
    const qd_real maxVal = static_cast<double> (std::numeric_limits<float>::max ());

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<float, qd_real>::safeConvert: "
      "Input qd_real t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to float.");

    // This should also test the double->float conversion for overflow.
    return asSafe<float> (to_double (t));
  }
};

/** \brief Convert qd_real to int. */
template <>
class ValueTypeConversionTraits<int, qd_real> {
public:
  inline static int convert (const qd_real t) {
    return to_int (t);
  }
  static int safeConvert (const qd_real t) {
    // qd_real has a constructor for int.
    const qd_real minVal = std::numeric_limits<int>::min ();
    const qd_real maxVal = std::numeric_limits<int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<int, qd_real>::safeConvert: "
      "Input qd_real t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to int.");
    return to_int (t);
  }
};

/** \brief Convert qd_real to dd_real. */
template <>
class ValueTypeConversionTraits<dd_real, qd_real> {
public:
  inline static dd_real convert (const qd_real t) {
    return to_dd_real(t);
  }
  static dd_real safeConvert (const qd_real t) {
    // std::numeric_limits<dd_real>::min() gives the minimum
    // _positive_ (normalized? not sure what this means for dd_real --
    // mfh 14 Nov 2012) value of type dd_real.  dd_real values are
    // built from two IEEE 754 doubles.  This means they can change
    // sign just by flipping the sign bit, so the "most negative"
    // finite dd_real is just the negative of the "most positive"
    // finite dd_real.
    //
    // qd_real has a constructor for dd_real.
    const qd_real minVal = -std::numeric_limits<dd_real>::max ();
    const qd_real maxVal = std::numeric_limits<dd_real>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<dd_real, qd_real>::safeConvert: "
      "Input qd_real t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to dd_real.");

    return to_dd_real (t);
  }
};

/** \brief Convert dd_real to double. */
template <>
class ValueTypeConversionTraits<double, dd_real> {
public:
  inline static double convert (const dd_real t) {
    return to_double (t);
  }
  static double safeConvert (const dd_real t) {
    // std::numeric_limits<double>::min() gives the minimum _positive_
    // normalized value of type double.  IEEE 754 floating-point
    // values can change sign just by flipping the sign bit, so the
    // "most negative" finite double is just the negative of the "most
    // positive" finite double.
    //
    // qd_real has a constructor for double.
    const dd_real minVal = -std::numeric_limits<double>::max ();
    const dd_real maxVal = std::numeric_limits<double>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<double, dd_real>::safeConvert: "
      "Input dd_real t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to double.");

    return to_double (t);
  }
};

/** \brief Convert dd_real to float. */
template <>
class ValueTypeConversionTraits<float, dd_real> {
public:
  inline static float convert (const dd_real t) {
    // This also checks for double->float overflow in a debug build.
    return as<float> (to_double (t));
  }
  static float safeConvert (const dd_real t) {
    // std::numeric_limits<float>::min() gives the minimum _positive_
    // normalized value of type float.  IEEE 754 floating-point
    // values can change sign just by flipping the sign bit, so the
    // "most negative" finite float is just the negative of the "most
    // positive" finite float.
    //
    // dd_real has a constructor for double but not for float,
    // so we cast to double first.
    const dd_real minVal = static_cast<double> (-std::numeric_limits<float>::max ());
    const dd_real maxVal = static_cast<double> (std::numeric_limits<float>::max ());

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<float, dd_real>::safeConvert: "
      "Input dd_real t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to float.");

    // This also checks for double->float overflow.
    return as<float> (to_double (t));
  }
};

/** \brief Convert dd_real to int. */
template <>
class ValueTypeConversionTraits<int, dd_real> {
public:
  inline static int convert (const dd_real t) {
    return to_int (t);
  }
  static int safeConvert (const dd_real t) {
    // dd_real has a constructor for int.
    const dd_real minVal = std::numeric_limits<int>::min ();
    const dd_real maxVal = std::numeric_limits<int>::max ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      t < minVal || t > maxVal,
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<int, dd_real>::safeConvert: "
      "Input dd_real t = " << t << " is out of the valid range [" << minVal
      << ", " << maxVal << "] for conversion to int.");
    return to_int (t);
  }
};

/** \brief Convert long unsigned int to qd_real. */
template <>
class ValueTypeConversionTraits<qd_real, long unsigned int> {
public:
  inline static qd_real convert( const long unsigned int t ) {
    // FIXME (mfh 14 Nov 2012): qd_real unfortunately lacks a
    // constructor or conversion function for conversions from
    // built-in integer types other than int.  However, it does allow
    // reading in values from a string.  We could use this to convert
    // from any type to qd_real, by first writing the value to an
    // std::ostringstream, then creating a qd_real from the resulting
    // string.
    return ValueTypeConversionTraits<qd_real,int>::convert(ValueTypeConversionTraits<int,long unsigned int>::convert(t));
  }
  inline static qd_real safeConvert( const long unsigned int t )
    { return ValueTypeConversionTraits<qd_real,int>::safeConvert(ValueTypeConversionTraits<int,long unsigned int>::safeConvert(t)); }
};

/** \brief Convert long unsigned int to dd_real. */
template <>
class ValueTypeConversionTraits<dd_real, long unsigned int> {
public:
  inline static dd_real convert( const long unsigned int t ) {
    // FIXME (mfh 14 Nov 2012): dd_real unfortunately lacks a
    // constructor or conversion function for conversions from
    // built-in integer types other than int.  However, it does allow
    // reading in values from a string.  We could use this to convert
    // from any type to dd_real, by first writing the value to an
    // std::ostringstream, then creating a dd_real from the resulting
    // string.
    return ValueTypeConversionTraits<dd_real,int>::convert(ValueTypeConversionTraits<int,long unsigned int>::convert(t));
  }
  inline static dd_real safeConvert( const long unsigned int t )
    { return ValueTypeConversionTraits<dd_real,int>::safeConvert(ValueTypeConversionTraits<int,long unsigned int>::safeConvert(t)); }
};

#endif // HAVE_TEUCHOS_QD

// ToDo: Add more specializations as needed!


} // end namespace Teuchos


#endif // TEUCHOS_AS_HPP
