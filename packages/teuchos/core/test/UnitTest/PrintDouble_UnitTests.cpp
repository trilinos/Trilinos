// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_PrintDouble.hpp"

#include <cstdlib>
#include <sstream>
#include <cmath>

namespace {

void test_print_double(double v, const char* s, bool& success, std::ostream& out) {
  std::stringstream ss;
  Teuchos::print_double(ss, v);
  auto s2 = ss.str();
  TEST_EQUALITY_CONST( s2, s );
  if (std::isnan(v)) {
    // NaN != NaN, special case
    TEST_EQUALITY_CONST( true, std::isnan(std::atof(s)) );
  } else {
    TEST_EQUALITY_CONST( v, std::atof(s) );
  }
}

void test_print_double(const char* s, bool& success, std::ostream& out) {
  test_print_double(std::atof(s), s, success, out);
}

double form_double(bool negative, std::uint64_t mantissa, std::int32_t exponent) {
  std::uint64_t pun = 0;
  if (negative) pun |= std::uint64_t(1) << 63;
  auto biased_exponent = std::uint64_t(exponent) + 1075;
  pun |= biased_exponent << 52;
  pun |= mantissa % (std::uint64_t(1) << 52);
  double value;
  std::memcpy(&value, &pun, sizeof(value));
  return value;
}

TEUCHOS_UNIT_TEST( PrintDouble, Basic )
{
  test_print_double("1.", success, out);
  test_print_double("-1.", success, out);
  test_print_double("0.", success, out);
  test_print_double("4.", success, out);
  test_print_double("55.", success, out);
  test_print_double("3.14159", success, out);
  test_print_double("5404319552844595.", success, out);
  test_print_double("1e300", success, out);
  test_print_double("1e-300", success, out);
  test_print_double("1e2", success, out);
  test_print_double("10.", success, out);
  test_print_double(".003", success, out);
  /* in order to ensure this gets converted back to the right
     value even if people change their IEEE rounding mode,
     this long representation is needed */
  test_print_double("9.999999999999999e22", success, out);
  test_print_double("inf", success, out);
  test_print_double("-inf", success, out);
  test_print_double("nan", success, out);
  // explicitly form 2^60 to test rare if case in print_double
  test_print_double(form_double(0, 0, 8), "1152921504606847000.", success, out);
  test_print_double(form_double(0, 0, 0), "4503599627370496.", success, out);
  test_print_double(form_double(0, 4503599627370503, -30), "4194304.0000000065", success, out);
  test_print_double(form_double(0, 4503599627370518, -30), "4194304.0000000205", success, out);
  // tests from Steele and White's paper
  test_print_double("1.3", success, out);
  test_print_double(4.0/3.0, "1.3333333333333333", success, out);
  test_print_double(1.0/10.0, ".1", success, out);
}

} // namespace
