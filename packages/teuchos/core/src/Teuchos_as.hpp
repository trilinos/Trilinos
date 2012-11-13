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
 * \section Teuchos_as_Dev Developer documentation
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

/// \brief Convert an std::string to a long, without compiler warnings.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<long, std::string> {
public:
  //! Convert the given std::string to a long.
  static long safeConvert (const std::string& t) {
    // We call strtol() instead of using std::istringstream, because
    // we want more detailed information in case of failure to
    // convert.  I have no idea what operator>>(std::istream&, long&)
    // does if it encounters an integer too long to fit in long, for
    // example.
    char* endptr = NULL;
    // Keep the pointer, because std::string doesn't necessarily
    // guarantee that this is the same across calls to c_str(), does
    // it?  Or perhaps it does...
    const char* t_ptr = t.c_str ();
    // We preset errno to 0, to distinguish success or failure after
    // calling strtol.  Most implementations of the C standard library
    // written with threads in mind have errno be a macro that expands
    // to thread-local storage.  Thanks to the Linux documentation for
    // strtol ("man 3 strtol", Red Hat Enterprise Linux 5) for advice
    // with the following checks.
    errno = 0;
    const long val = strtol (t_ptr, &endptr, 10);

    TEUCHOS_TEST_FOR_EXCEPTION(
      errno == ERANGE && (val == LONG_MAX || val == LONG_MIN),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<long, std::string>::convert: "
      "The integer value in the given std::string \"" << t << "\" overflows long.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      errno != 0 && val == 0,
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<long, std::string>::convert: "
      "stdtol was unable to convert the given std::string \"" << t << "\" to a long.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      endptr == t_ptr, // See above discussion of c_str().
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<long, std::string>::convert: "
      "stdtol was unable to read any integer digits from the given std::string \"" << t << "\".");

    return val;
  }

  //! Convert the given std::string to a long.
  static long convert (const std::string& t) {
    return safeConvert (t);
  }
};


/// \brief Convert an std::string to a unsigned long, without compiler warnings.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<unsigned long, std::string> {
public:
  //! Convert the given std::string to a unsigned long.
  static unsigned long safeConvert (const std::string& t) {
    // We call strtoul() instead of using std::istringstream, because
    // we want more detailed information in case of failure to
    // convert.  I have no idea what operator>>(std::istream&,
    // unsigned long&) does if it encounters an integer too long to
    // fit in unsigned long, for example.
    char* endptr = NULL;
    // Keep the pointer, because std::string doesn't necessarily
    // guarantee that this is the same across calls to c_str(), does
    // it?  Or perhaps it does...
    const char* t_ptr = t.c_str ();
    // We preset errno to 0, to distinguish success or failure after
    // calling strtoul.  Most implementations of the C standard
    // library written with threads in mind have errno be a macro that
    // expands to thread-local storage.  Thanks to the Linux
    // documentation for strtol ("man 3 strtol", Red Hat Enterprise
    // Linux 5) for advice with the following checks.
    errno = 0;
    const unsigned long val = strtoul (t_ptr, &endptr, 10);

    TEUCHOS_TEST_FOR_EXCEPTION(
      errno == ERANGE && (val == ULONG_MAX || val == static_cast<unsigned long> (0)),
      std::range_error,
      "Teuchos::ValueTypeConversionTraits<unsigned long, std::string>::convert: "
      "The integer value in the given std::string \"" << t << "\" overflows unsigned long.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      errno != 0 && val == 0,
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<unsigned long, std::string>::convert: "
      "stdtol was unable to convert the given std::string \"" << t << "\" to an unsigned long.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      endptr == t_ptr, // See above discussion of c_str().
      std::invalid_argument,
      "Teuchos::ValueTypeConversionTraits<unsigned long, std::string>::convert: "
      "stdtol was unable to read any integer digits from the given std::string \"" << t << "\".");

    return val;
  }

  //! Convert the given std::string to an unsigned long.
  static unsigned long convert (const std::string& t) {
    return safeConvert (t);
  }
};


/// \brief Convert an std::string to an int, without compiler warnings.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<int, std::string> {
private:
  //! Convert the given std::string to an intermediate long.
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
  //! Convert the given std::string to an int.
  static int safeConvert (const std::string& t) {
    return asSafe<int> (safeConvertToLong (t));
  }

  //! Convert the given std::string to an int.
  static int convert (const std::string& t) {
    return as<int> (safeConvertToLong (t));
  }
};


/// \brief Convert an std::string to an unsigned int, without compiler warnings.
///
/// We assume the string stores a base-10 integer, if it stores an integer at all.
template<>
class ValueTypeConversionTraits<unsigned int, std::string> {
private:
  //! Convert the given std::string to an intermediate unsigned long.
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
  //! Convert the given std::string to an unsigned int.
  static unsigned int safeConvert (const std::string& t) {
    return asSafe<unsigned int> (safeConvertToUnsignedLong (t));
  }

  //! Convert the given std::string to an unsigned int.
  static unsigned int convert (const std::string& t) {
    return as<unsigned int> (safeConvertToUnsignedLong (t));
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
    const int minInt = std::numeric_limits<int>::min ();
    const int maxInt = std::numeric_limits<int>::max ();

    // It's legal for sizeof (unsigned int) == 8, but rare.  In that
    // rare case, float certainly won't overflow unsigned int, but we
    // wouldn't want to do the check below, because it casts maxInt to
    // float.
    if (sizeof (unsigned int) <= sizeof (float)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < static_cast<float> (minInt) || t > static_cast<float> (maxInt),
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<int, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range ["
        << minInt << ", " << maxInt << "] for conversion to int.");
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
    const unsigned int minInt = 0; // Had better be, since it's unsigned.
    const unsigned int maxInt = std::numeric_limits<unsigned int>::max ();

    // It's legal for sizeof (unsigned int) == 8, but rare.  In that
    // rare case, float certainly won't overflow unsigned int, but we
    // wouldn't want to do the check below, because it casts maxInt to
    // float.
    if (sizeof (unsigned int) <= sizeof (float)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        t < static_cast<float> (minInt) || t > static_cast<float> (maxInt),
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<unsigned int, float>::safeConvert: "
        "Input float t = " << t << " is out of the valid range [" << minInt
        << ", " << maxInt << "] for conversion to unsigned int.");
    }
    return static_cast<unsigned int> (t);
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
  /// \warning <tt>unsigned long</tt> integer values may overflow
  ///   <tt>int</tt>, depending on your platform.  You should use
  ///   safeConvert() if you aren't sure that the given <tt>unsigned
  ///   long</tt> value fits in an <tt>int</tt>.
  static int convert (const unsigned long t) {
    // Implicit conversion from unsigned long to int may cause
    // compiler warnings, but static_cast does not.
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

    // The C++ standard promises that sizeof (int) <= sizeof (unsigned long).
    if (sizeof (int) == sizeof (unsigned long)) {
      // The two types have the same number of bits.  Thus,
      // two's-complement arithmetic means that if casting from
      // unsigned long to int results in a negative number, it
      // overflowed.  Otherwise, it didn't overflow (same number of
      // bits).
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<int> (t) < static_cast<int> (0),
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<int, unsigned long>::safeConvert: "
        "Input unsigned long t = " << t << " is out of the valid range ["
        << minInt << ", " << maxInt << "] for conversion to int.");
    }
    else { // sizeof (int) < sizeof (unsigned long).
      // t is unsigned, so it is >= 0 by definition.
      // Casting from int to unsigned long won't overflow in this case.
      TEUCHOS_TEST_FOR_EXCEPTION(
        t > static_cast<unsigned long> (maxInt),
        std::range_error,
        "Teuchos::ValueTypeConversionTraits<int, unsigned long>::safeConvert: "
        "Input unsigned long t = " << t << " is out of the valid range ["
        << minInt << ", " << maxInt << "] for conversion to int.");
    }

    // Implicit conversion from unsigned long to int may cause
    // compiler warnings, but static_cast does not.
    return static_cast<int> (t);
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
    const float minFloat = std::numeric_limits<float>::min ();
    const float maxFloat = std::numeric_limits<float>::max ();

    // t >= 0 by definition, because it is unsigned.
    TEUCHOS_TEST_FOR_EXCEPTION(
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
