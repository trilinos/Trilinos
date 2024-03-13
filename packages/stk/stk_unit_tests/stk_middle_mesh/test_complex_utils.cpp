#include <sstream>
#include "gtest/gtest.h"
#include "stk_middle_mesh/complex_utils.hpp"

namespace {

using Complex = std::complex<double>;

void expect_near(const Complex& a, const Complex& b, double tol)
{
  EXPECT_NEAR(a.real(), b.real(), tol);
  EXPECT_NEAR(a.imag(), b.imag(), tol);
}

}

TEST(ComplexUtils, Addition)
{
  Complex a(1, 2);
  expect_near(a + 2, Complex(3, 2), 1e-13);
  expect_near(2 + a, Complex(3, 2), 1e-13);
}

TEST(ComplexUtils, Subtraction)
{
  Complex a(1, 2);
  expect_near(a - 2, Complex(-1, 2), 1e-13);
  expect_near(2 - a, Complex( 1, -2), 1e-13);
}

TEST(ComplexUtils, Multiplication)
{
  Complex a(1, 2);
  expect_near(a * 2, Complex(2, 4), 1e-13);
  expect_near(2 * a, Complex(2, 4), 1e-13);
}

TEST(ComplexUtils, Division)
{
  Complex a(1, 2);
  expect_near(a / 2, Complex(0.5, 1), 1e-13);
  expect_near(2 / a, Complex(2.0/5, -4.0/5), 1e-13);
}


// This stuff is provided by the standard library, the tests are here just to
// verify it really works
TEST(ComplexUtils, Power)
{
  Complex a(1, 2);
  expect_near(std::pow(a, 2), Complex(-3, 4), 1e-13);
  expect_near(std::pow(a, 2.0), Complex(-3, 4), 1e-13);
}

TEST(ComplexUtils, OutputOperator)
{
  Complex a(1, 2);
  std::stringstream ss;
  ss << a;
  EXPECT_EQ(ss.str(), "(1,2)");


}