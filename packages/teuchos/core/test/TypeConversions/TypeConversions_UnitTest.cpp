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

#include "Teuchos_as.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <limits>

// Putting the unit tests in an anonymous namespace avoids name collisions.
namespace {

//
// Hack to work around Bug 5757 (unit test macros that instantiate
// templated unit tests can't handle spaces in the name of the
// template parameter).
//
typedef unsigned int unsigned_int_type;
typedef unsigned long unsigned_long_type;

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long long_long_type;
typedef unsigned long long unsigned_long_long_type;
#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// Templated test for overflow for conversion from double to a
// built-in signed integer type, if that should actually overflow
// (depends on sizeof(SignedIntType).
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, doubleToSignedIntTypeOverflow, SignedIntType )
{
  using Teuchos::asSafe;

  // std::numeric_limits<double>::min() gives the minimum _positive_
  // normalized value of type double.  IEEE 754 floating-point values
  // can change sign just by flipping the sign bit, so the "most
  // negative" finite double is just the negative of the "most
  // positive" finite double.
  const double minDouble = -std::numeric_limits<double>::max ();
  const double maxDouble = std::numeric_limits<double>::max ();

  SignedIntType val = 0;
  if (sizeof (SignedIntType) <= sizeof (double)) {
    TEST_THROW(val = asSafe<SignedIntType> (minDouble), std::range_error);
    TEST_THROW(val = asSafe<SignedIntType> (maxDouble), std::range_error);
  }
  (void) val; // Silence compiler errors.
}

//
// Templated test for overflow for conversion from double to a
// built-in unsigned integer type, if that should actually overflow
// (depends on sizeof(UnsignedIntType).
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, doubleToUnsignedIntTypeOverflow, UnsignedIntType )
{
  using Teuchos::asSafe;

  // std::numeric_limits<double>::min() gives the minimum _positive_
  // normalized value of type double.  IEEE 754 floating-point values
  // can change sign just by flipping the sign bit, so the "most
  // negative" finite double is just the negative of the "most
  // positive" finite double.
  const double minDouble = -std::numeric_limits<double>::max ();
  const double maxDouble = std::numeric_limits<double>::max ();

  UnsignedIntType val = 0;
  if (sizeof (UnsignedIntType) <= sizeof (double)) {
    TEST_THROW(val = asSafe<UnsignedIntType> (minDouble), std::range_error);
    TEST_THROW(val = asSafe<UnsignedIntType> (maxDouble), std::range_error);
    (void) val; // Silence compiler errors.
  }
  else { // Only conversions from negative values should throw.
    TEST_THROW(val = asSafe<UnsignedIntType> (minDouble), std::range_error);
    TEST_NOTHROW(val = asSafe<UnsignedIntType> (maxDouble));
    TEST_EQUALITY_CONST(val, static_cast<UnsignedIntType> (maxDouble));
  }

  // Conversion from any negative value should throw.
  const double minusOne = -1;
  TEST_THROW(val = asSafe<UnsignedIntType> (minusOne), std::range_error);
  (void) val; // Silence compiler errors.
}

//
// Instantiations of templated tests for conversions between various
// built-in integer types and double.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, doubleToSignedIntTypeOverflow, int )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, doubleToSignedIntTypeOverflow, long )

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, doubleToSignedIntTypeOverflow, long_long_type )
#endif // HAVE_TEUCHOS_LONG_LONG_INT

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, doubleToUnsignedIntTypeOverflow, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, doubleToUnsignedIntTypeOverflow, unsigned_long_type )

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, doubleToUnsignedIntTypeOverflow, unsigned_long_long_type )
#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// Templated test for conversions between possibly different built-in
// integer types.  The C++ standard guarantees the following:
// - <tt>sizeof (char) == 1</tt>
// - <tt>sizeof (char) <= sizeof (short) <= sizeof (int) <= sizeof (long)</tt>
//
// C99 also guarantees <tt>sizeof (long) <= sizeof (long long)</tt>.
//
// This means that any value that fits in a <tt>signed char</tt> must
// also fit in any other built-in integer type.  (The standard does
// not promise whether <tt>char</tt> is signed or unsigned.)
//
// We test both as() and asSafe to ensure correct behavior of both.
//

// Test for conversion between two built-in integer types FirstIntType
// and SecondIntType.  The test uses a positive number that must fit
// in both and must not overflow.  The test covers both as() and
// asSafe().
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( as, positiveFirstIntToSecondInt, FirstIntType, SecondIntType )
{
  using Teuchos::as;
  using Teuchos::asSafe;

  std::ostringstream os;
  const FirstIntType origVal = 42;
  const SecondIntType origValSecond = 42;

  SecondIntType asVal = 0, asSafeVal = 0;
  TEST_NOTHROW(asVal = as<SecondIntType> (origVal));
  TEST_NOTHROW(asSafeVal = asSafe<SecondIntType> (origVal));

  TEST_EQUALITY_CONST(asVal, static_cast<SecondIntType> (origValSecond));
  TEST_EQUALITY_CONST(asSafeVal, static_cast<SecondIntType> (origValSecond));
  TEST_EQUALITY_CONST(asVal, asSafeVal);

  FirstIntType backVal = 0, backSafeVal = 0;
  TEST_NOTHROW(backVal = as<FirstIntType> (asVal));
  TEST_NOTHROW(backSafeVal = asSafe<FirstIntType> (asSafeVal));

  TEST_EQUALITY_CONST(backVal, origVal);
  TEST_EQUALITY_CONST(backSafeVal, origVal);
  TEST_EQUALITY_CONST(backVal, backSafeVal);
}

// Test for conversion between two built-in integer types
// SignedIntType and UnsignedIntType.  The two types must have the
// same number of bits.  The test starts with a negative number that
// should trigger asSafe() to throw std::range_error.  as() will only
// throw std::range_error in a debug build, so we don't test it here.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( asSafe, negativeSignedIntToUnsignedInt, SignedIntType, UnsignedIntType )
{
  using Teuchos::asSafe;

  // Ensure that the two types have the same number of bits.
  TEUCHOS_TEST_FOR_EXCEPTION(
    sizeof (SignedIntType) != sizeof (UnsignedIntType),
    std::logic_error,
    "Unit test Teuchos,asSafe,negativeSignedIntToUnsignedInt requires that the "
    "two template parameters SignedIntType and UnsignedIntType have the same "
    "number of bits.");

  std::ostringstream os;
  const SignedIntType origVal = -1;

  UnsignedIntType asSafeVal = 0;
  // Casts from negative signed values to unsigned values should
  // throw, because they are not within range [0, maxUnsignedVal] of
  // the target type.
  TEST_THROW(asSafeVal = asSafe<UnsignedIntType> (origVal), std::range_error);
  (void) asSafeVal; // Silence compiler warning.

  // Casts from large unsigned values to negative signed values should
  // throw, because they change positivity of the result.
  UnsignedIntType negVal = static_cast<UnsignedIntType> (origVal);
  SignedIntType backSafeVal = 0;
  TEST_THROW(backSafeVal = asSafe<SignedIntType> (negVal), std::range_error);
  (void) backSafeVal; // Silence compiler warning.
}

// Test for conversion between two built-in signed integer types
// FirstSignedIntType and SecondSignedIntType.  The test uses a
// negative number that should not overflow in either case.  It tests
// both as() and asSafe().
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( as, negativeSignedIntToSignedInt, FirstSignedIntType, SecondSignedIntType )
{
  using Teuchos::as;
  using Teuchos::asSafe;

  std::ostringstream os;
  const FirstSignedIntType origVal = -42;

  // Ensure that the two types are both signed.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! std::numeric_limits<FirstSignedIntType>::is_signed ||
    ! std::numeric_limits<SecondSignedIntType>::is_signed,
    std::logic_error,
    "Unit test Teuchos,as,negativeSignedIntToSignedInt requires that the "
    "two template parameters FirstSignedIntType and SecondSignedIntType "
    "both be signed built-in integer types.");

  // Test cast from FirstSignedIntType to SecondSignedIntType.
  // The casts should not throw in either a debug or a release build.
  SecondSignedIntType asVal = 0, asSafeVal = 0;
  TEST_NOTHROW(asVal = as<SecondSignedIntType> (origVal));
  TEST_NOTHROW(asSafeVal = asSafe<SecondSignedIntType> (origVal));
  TEST_EQUALITY_CONST(asVal, static_cast<SecondSignedIntType> (origVal));
  TEST_EQUALITY_CONST(asSafeVal, static_cast<SecondSignedIntType> (origVal));
  TEST_EQUALITY_CONST(asVal, asSafeVal);

  FirstSignedIntType backVal = 0, backSafeVal = 0;
  TEST_NOTHROW(backVal = as<FirstSignedIntType> (origVal));
  TEST_NOTHROW(backSafeVal = asSafe<FirstSignedIntType> (origVal));
  TEST_EQUALITY_CONST(backVal, origVal);
  TEST_EQUALITY_CONST(backSafeVal, origVal);
  TEST_EQUALITY_CONST(backVal, backSafeVal);
}

//
// Instantiations of templated tests for conversions between two
// possibly different built-in integer types.
//

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, int )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, long )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, unsigned_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, unsigned_int_type )

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, long_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, long_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, long )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, long_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, int )
#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// Instantiations of templated tests for conversions from signed to
// unsigned built-in integer types.  The two types must have the same
// number of bits.
//

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, negativeSignedIntToUnsignedInt, int, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, negativeSignedIntToUnsignedInt, long, unsigned_long_type )

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, negativeSignedIntToUnsignedInt, long_long_type, unsigned_long_long_type )
#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// Instantiations of templated tests for conversions between two
// possibly different signed integer types, for a negative value that
// should not overflow.
//

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, int, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, int, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long, long )

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, int, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long_long_type, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long_long_type, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long_long_type, long_long_type )
#endif // HAVE_TEUCHOS_LONG_LONG_INT

//
// Tests for conversions from std::string to built-in integer types.
//

// Test for overflow when converting an std::string (containing a
// positive integer too large to fit in int) to int.
TEUCHOS_UNIT_TEST( asSafe, stringToIntPositiveOverflow ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  const int maxInt = std::numeric_limits<int>::max ();

  // Write a long int to the given string that is guaranteed to
  // overflow int, if long is strictly bigger than int.
  if (sizeof (int) < sizeof (long)) {
    const long maxIntPlusOne = static_cast<long> (maxInt) + static_cast<long> (1);
    os << maxIntPlusOne;

    // Now attempt to convert the string to int.  The conversion
    // should fail, but leave the string unaffected.
    int intVal = 0;
    TEST_THROW(intVal = asSafe<int> (os.str ()), std::range_error);
    (void) intVal; // Silence compiler warning.

    // Since the string is unaffected, conversion to long should work
    // just fine, and return the correct result.
    long longVal = 0;
    TEST_NOTHROW(longVal = asSafe<long> (os.str ()));
    TEST_EQUALITY_CONST(longVal, maxIntPlusOne);
  }
  else { // C++ standard requires then that sizeof(int) == sizeof(long).
    os << maxInt;
    // Converting the string to int should not throw and should return
    // the correct result.
    int intVal = 0;
    TEST_NOTHROW(intVal = asSafe<int> (os.str ()));
    TEST_EQUALITY_CONST(intVal, maxInt);
  }
}

// Test for overflow when converting an std::string (containing a
// negative integer too negative to fit in int) to int.
TEUCHOS_UNIT_TEST( asSafe, stringToIntNegativeOverflow ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  const int minInt = std::numeric_limits<int>::min ();

  // Write a long int to the given string that is guaranteed to
  // overflow int, if long is strictly bigger than int.
  if (sizeof (int) < sizeof (long)) {
    const long minIntMinusOne = static_cast<long> (minInt) - static_cast<long> (1);
    os << minIntMinusOne;

    // Now attempt to convert the string to int.  The conversion
    // should fail, but leave the string unaffected.
    int intVal = 0;
    TEST_THROW(intVal = asSafe<int> (os.str ()), std::range_error);
    (void) intVal; // Silence compiler warning

    // Since the string is unaffected, conversion to long should work
    // just fine, and return the correct result.
    long longVal = 0;
    TEST_NOTHROW(longVal = asSafe<long> (os.str ()));
    TEST_EQUALITY_CONST(longVal, minIntMinusOne);
  }
  else { // C++ standard requires then that sizeof(int) == sizeof(long).
    os << minInt;
    // Converting the string to int should not throw and should return
    // the correct result.
    int intVal = 0;
    TEST_NOTHROW(intVal = asSafe<int> (os.str ()));
    TEST_EQUALITY_CONST(intVal, minInt);
  }
}

// Unit test for conversion from std::string (containing a positive
// integer) to built-in integer types (may be signed or unsigned).
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, stringToIntegerPositive, IntegerType ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  os << static_cast<IntegerType> (42);
  IntegerType val = 0;
  TEST_NOTHROW(val = asSafe<IntegerType> (os.str ()));
  TEST_EQUALITY_CONST(val, static_cast<IntegerType> (42));
}

// Unit test for conversion from std::string (containing a negative
// integer) to built-in integer types (must be signed).
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, stringToIntegerNegative, IntegerType ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  os << static_cast<IntegerType> (-42);
  IntegerType val = 0;
  TEST_NOTHROW(val = asSafe<IntegerType> (os.str ()));
  TEST_EQUALITY_CONST(val, static_cast<IntegerType> (-42));
}

// Unit test for conversion from std::string (NOT containing an
// integer) to built-in integer types (may be signed or unsigned).
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, stringToIntegerShouldThrow, IntegerType ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  os << "This string definitely does not contain an integer.";
  IntegerType val = 0;
  TEST_THROW(val = asSafe<IntegerType> (os.str ()), std::invalid_argument);
  (void) val; // Silence compiler warning
}

// Macros to instantiate templated unit tests for conversion from
// std::string to built-in integer types.  AnyIntegerType may be
// signed or unsigned; SignedIntegerType must be signed.

#define UNIT_TEST_GROUP_ANY_INTEGER( AnyIntegerType ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, stringToIntegerPositive, AnyIntegerType ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, stringToIntegerShouldThrow, AnyIntegerType )

#define UNIT_TEST_GROUP_SIGNED_INTEGER( SignedIntegerType ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, stringToIntegerNegative, SignedIntegerType )

//
// Instantiations of templated unit tests for conversion from
// std::string to built-in integer types.
//

UNIT_TEST_GROUP_ANY_INTEGER( int )
UNIT_TEST_GROUP_SIGNED_INTEGER( int )
UNIT_TEST_GROUP_ANY_INTEGER( long )
UNIT_TEST_GROUP_SIGNED_INTEGER( long )

//UNIT_TEST_GROUP_ANY_INTEGER( unsigned int )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_int_type )
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned long )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_long_type )

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
//UNIT_TEST_GROUP_ANY_INTEGER( long long )
UNIT_TEST_GROUP_ANY_INTEGER( long_long_type )
//UNIT_TEST_GROUP_SIGNED_INTEGER( long long )
UNIT_TEST_GROUP_SIGNED_INTEGER( long_long_type )
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned long long )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_long_long_type )
#endif // HAVE_TEUCHOS_LONG_LONG_INT

} // namespace (anonymous)



