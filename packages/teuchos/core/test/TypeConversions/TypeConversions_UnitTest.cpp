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

// Like TEST_NOTHROW, but showing the exception message if the code does throw.
// The exception must be a subclass of std::exception.
#define TEST_NOTHROW_WITH_MESSAGE( code ) \
  try { \
    (out) << "Test that code {"#code";} does not throw : "; \
    code; \
    (out) << "passes\n"; \
  } \
  catch (std::exception& theException) {      \
    (success) = false; \
    out << "failed\n"; \
    out << "\nException message for unexpected exception:\n\n"; \
    { \
      Teuchos::OSTab l_tab(out); \
      out << theException.what() << "\n\n"; \
    } \
  }

// Putting the unit tests in an anonymous namespace avoids name collisions.
namespace {

//
// Hack to work around Bug 5757 (unit test macros that instantiate
// templated unit tests can't handle spaces in the name of the
// template parameter).
//
typedef unsigned short unsigned_short_type;
typedef unsigned int unsigned_int_type;
typedef unsigned long unsigned_long_type;

typedef long long long_long_type;
typedef unsigned long long unsigned_long_long_type;


template <class T>
std::string valToString(const T &val)
{
  const int precision = std::numeric_limits<T>::digits10 + 10; // Be super safe!
  //std::cout << "precision = " << precision << "\n";
  std::ostringstream os;
  os.precision(precision);
  os << val;
  return os.str();
}


//
// Tests for conversions between built-in floating-point types.
//
TEUCHOS_UNIT_TEST( asSafe, realToReal ) {
  using Teuchos::as;
  using Teuchos::asSafe;

  const float minF = -std::numeric_limits<float>::max ();
  const float minusOneF = -1;
  const float maxF = std::numeric_limits<float>::max ();

  const double minD = -std::numeric_limits<double>::max ();
  const double minusOneD = -1;
  const double maxD = std::numeric_limits<double>::max ();

  float valF = 0;
  double valD = 0;
  long double valLD = 0;

  //
  // Test float -> float conversions.
  //
  TEST_NOTHROW(valF = asSafe<float> (minF));
  TEST_EQUALITY_CONST(valF, minF);
  TEST_NOTHROW(valF = as<float> (minF));
  TEST_EQUALITY_CONST(valF, minF);
  TEST_NOTHROW(valF = asSafe<float> (maxF));
  TEST_EQUALITY_CONST(valF, maxF);
  TEST_NOTHROW(valF = as<float> (maxF));
  TEST_EQUALITY_CONST(valF, maxF);
  TEST_NOTHROW(valF = asSafe<float> (minusOneF));
  TEST_EQUALITY_CONST(valF, minusOneF);
  TEST_NOTHROW(valF = as<float> (minusOneF));
  TEST_EQUALITY_CONST(valF, minusOneF);

  //
  // Test double -> double conversions.
  //
  TEST_NOTHROW(valD = asSafe<double> (minD));
  TEST_EQUALITY_CONST(valD, minD);
  TEST_NOTHROW(valD = as<double> (minD));
  TEST_EQUALITY_CONST(valD, minD);
  TEST_NOTHROW(valD = asSafe<double> (maxD));
  TEST_EQUALITY_CONST(valD, maxD);
  TEST_NOTHROW(valD = as<double> (maxD));
  TEST_EQUALITY_CONST(valD, maxD);
  TEST_NOTHROW(valD = asSafe<double> (minusOneD));
  TEST_EQUALITY_CONST(valD, minusOneD);
  TEST_NOTHROW(valD = as<double> (minusOneD));
  TEST_EQUALITY_CONST(valD, minusOneD);

  //
  // Test double -> float conversions.
  //
  // mfh 25 Nov 2012: Disabled throwing std::range_error on overflow
  // for conversions between build-in floating-point types, in favor
  // of IEEE 754 overflow semantics.
  // TEST_THROW(valF = asSafe<float> (minD), std::range_error);
  // TEST_THROW(valF = asSafe<float> (maxD), std::range_error);

  TEST_NOTHROW(valF = asSafe<float> (minusOneF));
  TEST_EQUALITY_CONST(valF, minusOneF);
  TEST_NOTHROW(valF = as<float> (minusOneF));
  TEST_EQUALITY_CONST(valF, minusOneF);

  TEST_NOTHROW(valF = asSafe<float> (minusOneD));
  TEST_EQUALITY_CONST(valF, minusOneD);
  TEST_NOTHROW(valF = as<float> (minusOneD));
  TEST_EQUALITY_CONST(valF, minusOneD);

  //
  // Test float -> double conversions.
  //
  TEST_NOTHROW(valD = asSafe<double> (minF));
  TEST_EQUALITY_CONST(valD, minF);
  TEST_NOTHROW(valD = as<double> (minF));
  TEST_EQUALITY_CONST(valD, minF);

  TEST_NOTHROW(valD = asSafe<double> (maxF));
  TEST_EQUALITY_CONST(valD, maxF);
  TEST_NOTHROW(valD = as<double> (maxF));
  TEST_EQUALITY_CONST(valD, maxF);

  TEST_NOTHROW(valD = asSafe<double> (minusOneF));
  TEST_EQUALITY_CONST(valD, minusOneF);
  TEST_NOTHROW(valD = as<double> (minusOneF));
  TEST_EQUALITY_CONST(valD, minusOneF);

  // mfh 25 Nov 2012: C89 does not mandate that "long double"
  // implement the extended-precision 80-bit format of IEEE 754.  In
  // fact, Microsoft Visual Studio implements long double just as
  // double.  (This goes all the way back to Bill Gates' initial
  // discussions with the Intel x87 architects.)  Relaxing the sizeof
  // (long double) > sizeof (double) requirement will prevent test
  // failures such as the following (on Windows):
  //
  // http://testing.sandia.gov/cdash/testDetails.php?test=10628321&build=801972

  // // Make sure that long double is as long as the standard requires.
  // TEUCHOS_TEST_FOR_EXCEPTION(
  //   sizeof (long double) <= sizeof (double),
  //   std::logic_error,
  //   "Your system does not have an IEEE 754 - compliant implementation of long double.  "
  //   "The IEEE 754 standard requires that long double be longer than double.  "
  //   "In fact, it must use at least 80 bits. "
  //   "However, sizeof (long double) = " << sizeof (long double)
  //   << " < sizeof (double) = " << sizeof (double) << ".");

  const long double minLD = -std::numeric_limits<long double>::max ();
  const long double minusOneLD = -1;
  const long double maxLD = std::numeric_limits<long double>::max ();

  //
  // Test long double -> long double conversions.
  //
  TEST_NOTHROW(valLD = asSafe<long double> (minLD));
  TEST_EQUALITY_CONST(valLD, minLD);
  TEST_NOTHROW(valLD = as<long double> (minLD));
  TEST_EQUALITY_CONST(valLD, minLD);
  TEST_NOTHROW(valLD = asSafe<long double> (maxLD));
  TEST_EQUALITY_CONST(valLD, maxLD);
  TEST_NOTHROW(valLD = as<long double> (maxLD));
  TEST_EQUALITY_CONST(valLD, maxLD);
  TEST_NOTHROW(valLD = asSafe<long double> (minusOneLD));
  TEST_EQUALITY_CONST(valLD, minusOneLD);
  TEST_NOTHROW(valLD = as<long double> (minusOneLD));
  TEST_EQUALITY_CONST(valLD, minusOneLD);

  //
  // Test long double -> float conversions.
  //
  // mfh 25 Nov 2012: Disabled throwing std::range_error on overflow
  // for conversions between build-in floating-point types, in favor
  // of IEEE 754 overflow semantics.
  // TEST_THROW(valF = asSafe<float> (minLD), std::range_error);
  // TEST_THROW(valF = asSafe<float> (maxLD), std::range_error);
  TEST_NOTHROW(valF = asSafe<float> (minusOneLD));
  TEST_EQUALITY_CONST(valF, minusOneLD);
  TEST_NOTHROW(valF = as<float> (minusOneLD));
  TEST_EQUALITY_CONST(valF, minusOneLD);

  //
  // Test long double -> double conversions.
  //
  // mfh 25 Nov 2012: See above note on how long double is the same as
  // double on some systems.
  //
  // TEST_THROW(valD = asSafe<double> (minLD), std::range_error);
  // TEST_THROW(valD = asSafe<double> (maxLD), std::range_error);
  TEST_NOTHROW(valD = as<float> (minusOneLD));
  TEST_EQUALITY_CONST(valD, minusOneLD);

  //
  // Test float -> long double conversions.
  //
  TEST_NOTHROW(valLD = asSafe<long double> (minF));
  TEST_EQUALITY_CONST(valLD, minF);
  TEST_NOTHROW(valLD = as<long double> (minF));
  TEST_EQUALITY_CONST(valLD, minF);

  TEST_NOTHROW(valLD = asSafe<long double> (maxF));
  TEST_EQUALITY_CONST(valLD, maxF);
  TEST_NOTHROW(valLD = as<long double> (maxF));
  TEST_EQUALITY_CONST(valLD, maxF);

  TEST_NOTHROW(valLD = asSafe<long double> (minusOneF));
  TEST_EQUALITY_CONST(valLD, minusOneF);
  TEST_NOTHROW(valLD = as<long double> (minusOneF));
  TEST_EQUALITY_CONST(valLD, minusOneF);

  //
  // Test double -> long double conversions.
  //
  TEST_NOTHROW(valLD = asSafe<long double> (minD));
  TEST_EQUALITY_CONST(valLD, minD);
  TEST_NOTHROW(valLD = as<long double> (minD));
  TEST_EQUALITY_CONST(valLD, minD);

  TEST_NOTHROW(valLD = asSafe<long double> (maxD));
  TEST_EQUALITY_CONST(valLD, maxD);
  TEST_NOTHROW(valLD = as<long double> (maxD));
  TEST_EQUALITY_CONST(valLD, maxD);

  TEST_NOTHROW(valLD = asSafe<long double> (minusOneD));
  TEST_EQUALITY_CONST(valLD, minusOneD);
  TEST_NOTHROW(valLD = as<long double> (minusOneD));
  TEST_EQUALITY_CONST(valLD, minusOneD);
}


//
// Tests for conversions from std::string to built-in floating-point types.
//
TEUCHOS_UNIT_TEST( asSafe, stringToReal ) {
  using Teuchos::as;
  using Teuchos::asSafe;

  const float minF = -std::numeric_limits<float>::max ();
  const float minusOneF = -1;
  const float maxF = std::numeric_limits<float>::max ();

  out << "minF = " << minF << "\n";
  out << "maxF = " << maxF << "\n";

  const double minD = -std::numeric_limits<double>::max ();
  const double minusOneD = -1;
  const double maxD = std::numeric_limits<double>::max ();

  out << "minD = " << minD << "\n";
  out << "maxD = " << maxD << "\n";

  // mfh 26 Nov 2012: C89 does not mandate that "long double"
  // implement the extended-precision 80-bit format of IEEE 754.  In
  // fact, Microsoft Visual Studio implements long double just as
  // double.  (This goes all the way back to Bill Gates' initial
  // discussions with the Intel x87 architects.)  Relaxing the sizeof
  // (long double) > sizeof (double) requirement will prevent test
  // failures such as the following (on Windows):
  //
  // http://testing.sandia.gov/cdash/testDetails.php?test=10628321&build=801972
  // http://testing.sandia.gov/cdash/testDetails.php?test=10739503&build=810247

  // // Make sure that long double is as long as the standard requires.
  // TEUCHOS_TEST_FOR_EXCEPTION(
  //   sizeof (long double) <= sizeof (double),
  //   std::logic_error,
  //   "Your system does not have an IEEE 754 - compliant implementation of long double.  "
  //   "The IEEE 754 standard requires that long double be longer than double.  "
  //   "In fact, it must use at least 80 bits. "
  //   "However, sizeof (long double) = " << sizeof (long double)
  //   << " < sizeof (double) = " << sizeof (double) << ".");

  const long double minLD = -std::numeric_limits<long double>::max ();
  const long double minusOneLD = -1;
  const long double maxLD = std::numeric_limits<long double>::max ();

  out << "minLD = " << minLD << "\n";
  out << "maxLD = " << maxLD << "\n";

  float valF = 0;
  double valD = 0;
  long double valLD = 0;

  //
  out << "Testing string -> float conversions ...\n";
  //
  {
    std::ostringstream os;
    // mfh 27 Nov 2012: Write all 17 digits that the double deserves.
    // If you just write 9, it might round (as it does on some
    // platforms) to a value that can't be represented in float.
    // os.precision (9);
    os.precision (17);
    os << minF;
    TEST_NOTHROW_WITH_MESSAGE(valF = asSafe<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, minF);
    TEST_NOTHROW_WITH_MESSAGE(valF = as<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, minF);
  }
  {
    std::ostringstream os;
    os.precision (17);
    os << maxF;
    TEST_NOTHROW_WITH_MESSAGE(valF = asSafe<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, maxF);
    TEST_NOTHROW_WITH_MESSAGE(valF = as<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, maxF);
  }
  {
    std::ostringstream os;
    os.precision (17);
    os << minusOneF;
    TEST_NOTHROW_WITH_MESSAGE(valF = asSafe<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, minusOneF);
    TEST_NOTHROW_WITH_MESSAGE(valF = as<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, minusOneF);
  }
  // Write -1 as double, read as float; shouldn't throw.
  {
    std::ostringstream os;
    os.precision (17);
    os << minusOneD;
    TEST_NOTHROW_WITH_MESSAGE(valF = asSafe<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, minusOneF);
    TEST_NOTHROW_WITH_MESSAGE(valF = as<float> (os.str ()));
    TEST_EQUALITY_CONST(valF, minusOneF);
  }

  //
  out << "Testing string -> float conversions that should throw ...\n";
  //
  {
    std::ostringstream os;
    os.precision (9);
    os << minD;
    TEST_THROW(valF = asSafe<float> (os.str ()), std::range_error);
  }
  {
    std::ostringstream os;
    os.precision (9);
    os << maxD;
    TEST_THROW(valF = asSafe<float> (os.str ()), std::range_error);
  }

  //
  out << "Testing string -> double conversions ...\n";
  //
  {
    std::ostringstream os;
    os.precision (17);
    os << minD;
    TEST_NOTHROW(valD = asSafe<double> (os.str ()));
    TEST_EQUALITY_CONST(valD, minD);
    TEST_NOTHROW(valD = as<double> (os.str ()));
    TEST_EQUALITY_CONST(valD, minD);
  }
  {
    std::ostringstream os;
    os.precision (17);
    os << maxD;
    TEST_NOTHROW(valD = asSafe<double> (os.str ()));
    TEST_EQUALITY_CONST(valD, maxD);
    TEST_NOTHROW(valD = as<double> (os.str ()));
    TEST_EQUALITY_CONST(valD, maxD);
  }
  {
    TEST_NOTHROW(valD = asSafe<double> (valToString(minusOneD)));
    TEST_EQUALITY_CONST(valD, minusOneD);
    TEST_NOTHROW(valD = as<double> (valToString(minusOneD)));
    TEST_EQUALITY_CONST(valD, minusOneD);
  }

  //
  // Test string -> double conversions that should throw,
  // if sizeof(long double) > sizeof(double).
  //
  const int sizeof_long_double = sizeof(long double);
  const int sizeof_double = sizeof(double);
  const int max_exponent10_long_double = std::numeric_limits<long double>::max_exponent10;
  const int max_exponent10_double = std::numeric_limits<double>::max_exponent10;

  out << "sizeof_long_double = " << sizeof_long_double << "\n";
  out << "sizeof_double = " << sizeof_double << "\n";
  out << "max_exponent10_long_double = " << max_exponent10_long_double << "\n";
  out << "max_exponent10_double = " << max_exponent10_double << "\n";
 
  if (sizeof_long_double > sizeof_double
    && max_exponent10_long_double > max_exponent10_double)
  {
    out << "Testing converting from 'long double' to 'double' that does not fit ...\n";

    const long double maxTooBigLD_for_D = static_cast<long double>(maxD) * 10.0;
    const long double minTooBigLD_for_D = static_cast<long double>(minD) * 10.0;
    out << "maxTooBigLD_for_D = " << maxTooBigLD_for_D << "\n";
    out << "minTooBigLD_for_D = " << minTooBigLD_for_D << "\n";

    {
      TEST_THROW(valD = asSafe<double>(valToString(minTooBigLD_for_D)), std::range_error);
      out << "valD = " << valD << "\n";
    }
    {
      std::ostringstream os;
      os.precision (36);
      TEST_THROW(valD = asSafe<double>(valToString(maxTooBigLD_for_D)), std::range_error);
      out << "valD = " << valD << "\n";
    }

    // NOTE: The above test avoids using std::numeric_limits<long
    // double>::max() because with the CUDA compiler on shiller/hansen it
    // returns 'inf'.  That completely breaks the test since "inf" is a valid
    // value to read into a 'double'.
    //
    // This updated test takes std::numeric_limits<double>::max() * 10.0,
    // writes to a string and then tries to read that back in as a 'dobule'.
    // That cathces the bad conversion and thowns an exception as it should.
    // This will work on any system where std::numeric_limits<long
    // double>::max_expoent10 > std::numeric_limits<double>::max_expoent10
    // (and this test is only run in cases where that is true).  See Trilinos
    // GitHub issue #2407.

  }

  //
  out << "Testing string -> long double conversions ...\n";
  //
  {
    std::ostringstream os;
    os.precision (36);
    os << minLD;
    TEST_NOTHROW(valLD = asSafe<long double> (os.str ()));
    TEST_EQUALITY_CONST(valLD, minLD);
    TEST_NOTHROW(valLD = as<long double> (os.str ()));
    TEST_EQUALITY_CONST(valLD, minLD);
  }
  {
    std::ostringstream os;
    os.precision (36);
    os << maxLD;
    TEST_NOTHROW(valLD = asSafe<long double> (os.str ()));
    TEST_EQUALITY_CONST(valLD, maxLD);
    TEST_NOTHROW(valLD = as<long double> (os.str ()));
    TEST_EQUALITY_CONST(valLD, maxLD);
  }
  {
    std::ostringstream os;
    os.precision (36);
    os << minusOneLD;
    TEST_NOTHROW(valLD = asSafe<long double> (os.str ()));
    TEST_EQUALITY_CONST(valLD, minusOneLD);
    TEST_NOTHROW(valLD = as<long double> (os.str ()));
    TEST_EQUALITY_CONST(valLD, minusOneLD);
  }
}


//
// Templated test for overflow for conversion from a built-in
// real-valued (not complex) floating-point type (float or double) to
// a built-in signed integer type, if that should actually overflow
// (depends on sizeof(SignedIntType)).
//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( asSafe, realToSignedIntTypeOverflow, RealType, SignedIntType )
{
  using Teuchos::asSafe;

  // std::numeric_limits<RealType>::min() gives the minimum _positive_
  // normalized value of type RealType.  IEEE 754 floating-point
  // values can change sign just by flipping the sign bit, so the
  // "most negative" finite RealType is just the negative of the "most
  // positive" finite RealType.
  const RealType minVal = -std::numeric_limits<RealType>::max ();
  const RealType maxVal = std::numeric_limits<RealType>::max ();

  SignedIntType val = 0;
  if (sizeof (SignedIntType) < sizeof (RealType)) {
    TEST_THROW(val = asSafe<SignedIntType> (minVal), std::range_error);
    TEST_THROW(val = asSafe<SignedIntType> (maxVal), std::range_error);
  }
  (void) val; // Silence compiler errors.
}

//
// Templated test for overflow for conversion from a built-in
// real-valued (not complex) floating-point type (float or double) to
// a built-in unsigned integer type, if that should actually overflow
// (depends on sizeof(UnsignedIntType)).
//
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( asSafe, realToUnsignedIntTypeOverflow, RealType, UnsignedIntType )
{
  using Teuchos::asSafe;
  using Teuchos::TypeNameTraits;

  // std::numeric_limits<RealType>::min() gives the minimum _positive_
  // normalized value of type RealType.  IEEE 754 floating-point
  // values can change sign just by flipping the sign bit, so the
  // "most negative" finite RealType is just the negative of the "most
  // positive" finite RealType.
  const RealType minVal = -std::numeric_limits<RealType>::max ();
  const RealType maxVal = std::numeric_limits<RealType>::max ();
  const UnsignedIntType maxUnsignedIntVal =
    std::numeric_limits<UnsignedIntType>::max ();

  // mfh 15 Nov 2012: Set val to a marker value, so we can see if the
  // body of TEST_NOTHROW below actually did the assignment.
  UnsignedIntType val = 42;
  // Make sure that val starts off with a different value than what
  // its final value should be.
  TEUCHOS_TEST_FOR_EXCEPTION(
    val == static_cast<UnsignedIntType> (maxVal),
    std::logic_error,
    "Dear test author, please pick a different marker value.  "
    "Please report this bug to the Teuchos developers.");

  // Conversion from any negative value should throw.
  TEST_THROW(val = asSafe<UnsignedIntType> (minVal), std::range_error);
  const RealType minusOne = -1;
  TEST_THROW(val = asSafe<UnsignedIntType> (minusOne), std::range_error);

  // Only test overflow checks if overflow can actually take place.
  if (maxUnsignedIntVal < maxVal) {
    TEST_THROW(val = asSafe<UnsignedIntType> (maxVal), std::range_error);
    try {
      std::cerr << std::endl
		<< "*** RealType = " << TypeNameTraits<RealType>::name ()
		<< ", UnsignedIntType = " << TypeNameTraits<UnsignedIntType>::name ()
		<< ", maxVal = " << maxVal
                << ", maxUnsignedIntVal = " << maxUnsignedIntVal
		<< ", asSafe (maxVal) = " << asSafe<UnsignedIntType> (maxVal)
		<< std::endl;
    } catch (...) {
      std::cerr << "(asSafe threw an exception)" << std::endl;
    }
  }
  else { // Only conversions from negative values should throw.
    TEST_NOTHROW(val = asSafe<UnsignedIntType> (maxVal));
    TEST_EQUALITY_CONST(val, static_cast<UnsignedIntType> (maxVal));

#if 0
    TEUCHOS_TEST_FOR_EXCEPTION(
      val == 42,
      std::logic_error,
      "Hey, how come val == 42?  It should be something completely different.  "
      << std::endl
      << "FYI, static_cast<" << TypeNameTraits<UnsignedIntType>::name ()
      << "> (minVal) = " << static_cast<UnsignedIntType> (minVal)
      << " and "
      << std::endl
      << "static_cast<" << TypeNameTraits<UnsignedIntType>::name ()
      << "> (maxVal) = " << static_cast<UnsignedIntType> (maxVal)
      << ".  val should be equal to the latter."
      << std::endl
      << "As float: minVal = " << minVal << ", maxVal = " << maxVal << ".");
#endif // 0
  }

  (void) val; // Silence compiler errors.
}

//
// Instantiations of templated tests for conversions from double to
// various built-in integer types.
//

// mfh 19 Nov 2012: The tests that I disabled below in commit
// f99c0e446f5c8dc385d00b60878314d40a7b9fe2 appear to be working now.
// I am reenabling them tentatively.
//
// mfh 16 Nov 2012: The (asSafe, realToUnsignedIntTypeOverflow) test
// keeps failing for template parameter combinations (double,
// unsigned_long_type), (float, unsigned_int_type), and (float,
// unsigned_long_type).  It only fails on some platforms, not all, and
// I can't figure out why.  I'm disabling these tests for now until I
// get more time to deal with this.  For examples of test output, see:
//
// http://testing.sandia.gov/cdash/testDetails.php?test=10519753&build=793648
// http://testing.sandia.gov/cdash/testDetails.php?test=10519852&build=793655
// http://testing.sandia.gov/cdash/testDetails.php?test=10520495&build=793698
// http://testing.sandia.gov/cdash/testDetails.php?test=10523690&build=793963
// http://testing.sandia.gov/cdash/testDetails.php?test=10523763&build=793962
// http://testing.sandia.gov/cdash/testDetails.php?test=10530499&build=794533
// http://testing.sandia.gov/cdash/testDetails.php?test=10530585&build=794532
// http://testing.sandia.gov/cdash/testDetails.php?test=10535648&build=794860

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, double, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, double, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, double, long )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, double, long_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, double, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToUnsignedIntTypeOverflow, double, unsigned_int_type )
// mfh 16,19 Nov 2012: See note above on formerly disabled tests.
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToUnsignedIntTypeOverflow, double, unsigned_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToUnsignedIntTypeOverflow, double, unsigned_long_long_type )

//
// Instantiations of templated tests for conversions from float to
// various built-in integer types.
//

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, float, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, float, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, float, long )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, float, long_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToSignedIntTypeOverflow, float, unsigned_short_type )
// mfh 16,19 Nov 2012: See note above on formerly disabled tests.
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToUnsignedIntTypeOverflow, float, unsigned_int_type )
// mfh 16,19 Nov 2012: See note above on formerly disabled tests.
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToUnsignedIntTypeOverflow, float, unsigned_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, realToUnsignedIntTypeOverflow, float, unsigned_long_long_type )


//
// Templated test for conversions between possibly different built-in
// integer types.  The C++ standard guarantees the following:
// - <tt>sizeof (char) == 1</tt>
// - <tt>sizeof (char) <= sizeof (short) <= sizeof (int) <= sizeof (long)</tt>
//
// C99 actually guarantees minimum sizes of the various types, and
// that <tt>sizeof (long) <= sizeof (long long)</tt>.
//
// This means that any value that fits in a <tt>signed char</tt> must
// also fit in any other built-in integer type.  (The standard does
// not promise whether <tt>char</tt> is signed or unsigned.)  We've
// chosen test values accordingly.
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

//
// 1. Tests for types of the same size.
//
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, short, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_short_type, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, short, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_short_type, short )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, int )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, long )

//
// 2. Tests for types of possibly different sizes.
//
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, short, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_short_type, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, short, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_short_type, unsigned_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, int, unsigned_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_int_type, unsigned_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_type, unsigned_int_type )

//
// 3. "long long", "unsigned long long" tests
//
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, long_long_type )

//
// 4. Tests between "long long" or "unsigned long long", and some
//    other built-in integer type.
//
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

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, short, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_short_type, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, short, unsigned_long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_short_type, long_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, long_long_type, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, positiveFirstIntToSecondInt, unsigned_long_long_type, short )

//
// Instantiations of templated tests for conversions from signed to
// unsigned built-in integer types.  The two types must have the same
// number of bits.
//

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, negativeSignedIntToUnsignedInt, short, unsigned_short_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, negativeSignedIntToUnsignedInt, int, unsigned_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, negativeSignedIntToUnsignedInt, long, unsigned_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( asSafe, negativeSignedIntToUnsignedInt, long_long_type, unsigned_long_long_type )

//
// Instantiations of templated tests for conversions between two
// possibly different signed integer types, for a negative value that
// should not overflow.
//

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, short, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, short, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, short, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, int, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, int, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, int, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long, long )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, short, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, int, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long, long_long_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long_long_type, short )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long_long_type, int )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long_long_type, long )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( as, negativeSignedIntToSignedInt, long_long_type, long_long_type )

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

UNIT_TEST_GROUP_ANY_INTEGER( short )
UNIT_TEST_GROUP_SIGNED_INTEGER( short )
UNIT_TEST_GROUP_ANY_INTEGER( int )
UNIT_TEST_GROUP_SIGNED_INTEGER( int )
UNIT_TEST_GROUP_ANY_INTEGER( long )
UNIT_TEST_GROUP_SIGNED_INTEGER( long )

//UNIT_TEST_GROUP_ANY_INTEGER( unsigned short )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_short_type )
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned int )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_int_type )
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned long )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_long_type )

//UNIT_TEST_GROUP_ANY_INTEGER( long long )
UNIT_TEST_GROUP_ANY_INTEGER( long_long_type )
//UNIT_TEST_GROUP_SIGNED_INTEGER( long long )
UNIT_TEST_GROUP_SIGNED_INTEGER( long_long_type )
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned long long )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_long_long_type )

} // namespace (anonymous)



