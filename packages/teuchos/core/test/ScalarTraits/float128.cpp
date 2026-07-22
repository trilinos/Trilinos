// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Array.hpp"
#include "Teuchos_ScalarTraits.hpp" // operator<< and operator>> overloads
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace { // (anonymous)

using std::endl;

TEUCHOS_UNIT_TEST( Float128, OutputStreamOp )
{
  Teuchos::OSTab tab0 (out);

#ifdef HAVE_TEUCHOSCORE_QUADMATH
  out << "Test operator<< (std::ostream&, __float128)" << endl;

  __float128 x = 1.0;
  __float128 y = strtoflt128 ("1.111112222233333", NULL);
  // __float128 has enough digits to represent this exactly, but
  // double precision would round.
  __float128 z = strtoflt128 ("1.111112222233333444445555566666", NULL);

  // FIXME (mfh 04 Sep 2015) The results of printing could depend on
  // the locale.  This works fine for the default locale on my system.
  {
    std::ostringstream os;
    os << x;
    TEST_EQUALITY_CONST( os.str (), "1.000000000000000000000000000000e+00" );
  }
  {
    std::ostringstream os;
    os << y;
    TEST_EQUALITY_CONST( os.str (), "1.111112222233333000000000000000e+00" );
  }
  {
    std::ostringstream os;
    os << z;
    TEST_EQUALITY_CONST( os.str (), "1.111112222233333444445555566666e+00" );
  }

  // Test that operator<< (std::ostream&, __float128) works.  The
  // operator<< overload MUST be defined in the std namespace.  If
  // defined in the global namespace, the compiler will have trouble
  // with the TEST_COMPARE_ARRAYS expression below.

  out << "Test chaining operator<<:" << endl
      << "z = " << z << ", it really does." << endl;

  // Test op<<, but make the __float128 arguments entries of an array.
  // (Make sure that this compiles without "ambiguous overload"
  // errors.)
  __float128 a[3];
  a[0] = x;
  a[1] = y;
  a[2] = z;
  out << "Testing chaining operator<< with array entries:" << endl
      << "a[0] = " << a[0] << ", a[1] = " << a[1] << ", a[2] = " << a[2]
      << endl;

  Teuchos::Array<__float128> arrayOfFloat128 (1);
  Teuchos::Tuple<__float128, 1> tupleOfFloat128;
  arrayOfFloat128[0] = z;
  tupleOfFloat128[0] = z;
  Teuchos::ArrayView<__float128> arrayViewOfFloat128 = arrayOfFloat128 ();
  // mfh 04 Sep 2015: If operator<< (std::ostream&, __float128) is
  // defined in the global namespace, instead of in the std namespace,
  // the TEST_COMPARE_ARRAYS expression will fail to compile.  GCC
  // 5.2.0 complains about the following line of
  // Teuchos::compareArrays (in
  // teuchos/core/src/Teuchos_TestingHelpers.hpp):
  //
  // out << "\nError, "<<a1_name<<"["<<i<<"] = "<<a1[i]<<" == "
  //     << a2_name<<"["<<i<<"] = "<<a2[i]<<": failed!\n";
  //
  // where a1 is a Teuchos::ArrayView<__float128> and a2 is a
  // Teuchos::Tuple<__float128, 1>.  The compiler claims an ambiguous
  // overload.  It has something to do with compareArrays (and
  // analogous functions in that file) being templated on Array types,
  // as I don't see that when I directly imitate what that line of
  // code does (see below).

  TEST_COMPARE_ARRAYS( arrayViewOfFloat128, tupleOfFloat128 );

  std::string s1 ("Some string");
  out << "Hello there! \"" << s1 << "\" is a string that doesn't mean anything "
    "on its own, but just makes this line of code more like the line of code "
    "that doesn't compile.  arrayViewOfFloat128[0] = " << arrayViewOfFloat128[0]
      << " and tupleOfFloat128[0] = " << tupleOfFloat128[0] << "." << endl;

  // Test that operator>> (std::istream&, __float128&) works.
  {
    // Use enough digits in the example to catch "cheats" that
    // truncate to double.
    const std::string z_str ("1.111112222233333444445555566666");
    std::istringstream is (z_str);
    __float128 z_copy = 0.0;
    is >> z_copy;
    TEST_EQUALITY_CONST( z, z_copy );
  }

#else
  out << "This test only makes sense to run with libquadmath enabled." << endl;
#endif // HAVE_TEUCHOSCORE_QUADMATH
}

} // namespace (anonymous)
