//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_NumericalUtils.hpp"

namespace Tempus_Unit_Test {

using Tempus::approxEqual;
using Tempus::approxEqualAbsTol;
using Tempus::approxEqualScale;

static double PI  = M_PI;
static double eps = 1.0e-14;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NumericalUtils, approxEqual)
{
  double zero = 0.0;
  double one  = 1.0;

  std::vector<double> scales;
  scales.push_back(1.0);
  scales.push_back(100.0);
  scales.push_back(0.01);
  scales.push_back(1.0e-28);
  scales.push_back(1.0e+28);

  std::vector<double> numbers;
  numbers.push_back(one);
  numbers.push_back(PI);

  for (auto& s : scales) {
    out << "  *****************************" << std::endl;
    out << "  s, eps = " << s << ", " << eps << std::endl;

    TEST_COMPARE(approxEqual(zero * s, (zero - 10.0 * eps) * s), ==, false);
    TEST_COMPARE(approxEqual(zero * s, (zero - eps) * s), ==, false);
    TEST_COMPARE(approxEqual(zero * s, (zero)*s), ==, true);
    TEST_COMPARE(approxEqual(zero * s, (zero + eps) * s), ==, false);
    TEST_COMPARE(approxEqual(zero * s, (zero + 10.0 * eps) * s), ==, false);

    // Swap
    TEST_COMPARE(approxEqual((zero - 10.0 * eps) * s, zero * s), ==, false);
    TEST_COMPARE(approxEqual((zero - eps) * s, zero * s), ==, false);
    TEST_COMPARE(approxEqual((zero)*s, zero * s), ==, true);
    TEST_COMPARE(approxEqual((zero + eps) * s, zero * s), ==, false);
    TEST_COMPARE(approxEqual((zero + 10.0 * eps) * s, zero * s), ==, false);

    for (auto& n : numbers) {
      out << "  n = " << n << std::endl;

      TEST_COMPARE(approxEqual(n * s, (n - 10.0 * eps) * s), ==, false);
      TEST_COMPARE(approxEqual(n * s, (n - eps) * s), ==, true);
      TEST_COMPARE(approxEqual(n * s, (n)*s), ==, true);
      TEST_COMPARE(approxEqual(n * s, (n + eps) * s), ==, true);
      TEST_COMPARE(approxEqual(n * s, (n + 10.0 * eps) * s), ==, false);

      // Swap
      TEST_COMPARE(approxEqual((n - 10.0 * eps) * s, n * s), ==, false);
      TEST_COMPARE(approxEqual((n - eps) * s, n * s), ==, true);
      TEST_COMPARE(approxEqual((n)*s, n * s), ==, true);
      TEST_COMPARE(approxEqual((n + eps) * s, n * s), ==, true);
      TEST_COMPARE(approxEqual((n + 10.0 * eps) * s, n * s), ==, false);
    }
    out << "  *****************************" << std::endl;
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NumericalUtils, approxEqualAbsTol)
{
  double numTol = Tempus::numericalTol<double>();
  double zero   = 0.0;
  double one    = 1.0;

  std::vector<double> scales;
  scales.push_back(1.0);
  scales.push_back(100.0);
  scales.push_back(0.01);
  scales.push_back(1.0e-28);
  scales.push_back(1.0e+28);

  std::vector<double> numbers;
  numbers.push_back(zero);
  numbers.push_back(one);
  numbers.push_back(PI);

  for (auto& s : scales) {
    for (auto& n : numbers) {
      out << "  *****************************" << std::endl;
      out << "  n, s, eps = " << n << ", " << s << ", " << eps << std::endl;

      TEST_COMPARE(approxEqualAbsTol(n * s, (n - 10.0 * eps) * s, numTol * s),
                   ==, false);
      TEST_COMPARE(approxEqualAbsTol(n * s, (n - eps) * s, numTol * s), ==,
                   true);
      TEST_COMPARE(approxEqualAbsTol(n * s, (n)*s, numTol * s), ==, true);
      TEST_COMPARE(approxEqualAbsTol(n * s, (n + eps) * s, numTol * s), ==,
                   true);
      TEST_COMPARE(approxEqualAbsTol(n * s, (n + 10.0 * eps) * s, numTol * s),
                   ==, false);

      // Swap
      TEST_COMPARE(approxEqualAbsTol((n - 10.0 * eps) * s, n * s, numTol * s),
                   ==, false);
      TEST_COMPARE(approxEqualAbsTol((n - eps) * s, n * s, numTol * s), ==,
                   true);
      TEST_COMPARE(approxEqualAbsTol((n)*s, n * s, numTol * s), ==, true);
      TEST_COMPARE(approxEqualAbsTol((n + eps) * s, n * s, numTol * s), ==,
                   true);
      TEST_COMPARE(approxEqualAbsTol((n + 10.0 * eps) * s, n * s, numTol * s),
                   ==, false);

      out << "  *****************************" << std::endl;
    }
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NumericalUtils, approxEqualScale)
{
  double zero = 0.0;
  double one  = 1.0;

  std::vector<double> scales;
  scales.push_back(1.0);
  scales.push_back(100.0);
  scales.push_back(0.01);
  scales.push_back(1.0e-28);
  scales.push_back(1.0e+28);

  std::vector<double> numbers;
  numbers.push_back(zero);
  numbers.push_back(one);
  numbers.push_back(PI);

  for (auto& s : scales) {
    for (auto& n : numbers) {
      out << "  *****************************" << std::endl;
      out << "  n, s, eps = " << n << ", " << s << ", " << eps << std::endl;

      TEST_COMPARE(approxEqualScale(n * s, (n - 10.0 * eps) * s, s), ==, false);
      TEST_COMPARE(approxEqualScale(n * s, (n - eps) * s, s), ==, true);
      TEST_COMPARE(approxEqualScale(n * s, (n)*s, s), ==, true);
      TEST_COMPARE(approxEqualScale(n * s, (n + eps) * s, s), ==, true);
      TEST_COMPARE(approxEqualScale(n * s, (n + 10.0 * eps) * s, s), ==, false);

      // Swap
      TEST_COMPARE(approxEqualScale((n - 10.0 * eps) * s, n * s, s), ==, false);
      TEST_COMPARE(approxEqualScale((n - eps) * s, n * s, s), ==, true);
      TEST_COMPARE(approxEqualScale((n)*s, n * s, s), ==, true);
      TEST_COMPARE(approxEqualScale((n + eps) * s, n * s, s), ==, true);
      TEST_COMPARE(approxEqualScale((n + 10.0 * eps) * s, n * s, s), ==, false);

      out << "  *****************************" << std::endl;
    }
  }
}

}  // namespace Tempus_Unit_Test
