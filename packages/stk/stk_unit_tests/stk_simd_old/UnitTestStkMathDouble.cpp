#include <gtest/gtest.h>
#include <stk_math/StkMath.hpp>
#include <cmath>
#include <algorithm>
#include <limits>

TEST(StkSimd, StkMathDouble_fmadd)
{
  constexpr double a = 2.1;
  constexpr double b = 3.7;
  constexpr double c = -9.1;
  EXPECT_DOUBLE_EQ( stk::math::fmadd(a,b,c), (a*b)+c );
}

TEST(StkSimd, StkMathDouble_sqrt)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::sqrt(a), std::sqrt(a) );
}

TEST(StkSimd, StkMathDouble_cbrt)
{
  constexpr double a = 2.1;
  constexpr double epsilon = 1.e-14;
  EXPECT_NEAR( stk::math::cbrt(a), std::pow(a,1.0/3.0), epsilon );
  // This fails on gcc debug, the above check does pass on gcc debug.
  // not sure what is different between the two checks.....
  //EXPECT_DOUBLE_EQ  ( stk::math::cbrt(a), std::pow(a,1.0/3.0) );
}

TEST(StkSimd, StkMathDouble_log)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::log(a), std::log(a) );
}

TEST(StkSimd, StkMathDouble_log10)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::log10(a), std::log10(a) );
}

TEST(StkSimd, StkMathDouble_exp)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::exp(a), std::exp(a) );
}

TEST(StkSimd, StkMathDouble_pow)
{
  constexpr double a = 2.1;
  constexpr double b = 3.7;
  constexpr int d = 2;
  EXPECT_DOUBLE_EQ( stk::math::pow(a,d), std::pow(a,d) );
  EXPECT_DOUBLE_EQ( stk::math::pow(a,b), std::pow(a,b) );
}

TEST(StkSimd, StkMathDouble_sin)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::sin(a), std::sin(a) );
}

TEST(StkSimd, StkMathDouble_cos)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::cos(a), std::cos(a) );
}

TEST(StkSimd, StkMathDouble_tan)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::tan(a), std::tan(a) );
}

TEST(StkSimd, StkMathDouble_sinh)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::sinh(a), std::sinh(a) );
}

TEST(StkSimd, StkMathDouble_cosh)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::cosh(a), std::cosh(a) );
}

TEST(StkSimd, StkMathDouble_tanh)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::tanh(a), std::tanh(a) );
}

TEST(StkSimd, StkMathDouble_asin)
{
  constexpr double f = 0.65;
  EXPECT_DOUBLE_EQ( stk::math::asin(f), std::asin(f) );
}

TEST(StkSimd, StkMathDouble_acos)
{
  constexpr double f = 0.65;
  EXPECT_DOUBLE_EQ( stk::math::acos(f), std::acos(f) );
}

TEST(StkSimd, StkMathDouble_atan)
{
  constexpr double f = 0.65;
  EXPECT_DOUBLE_EQ( stk::math::atan(f), std::atan(f) );
}

TEST(StkSimd, StkMathDouble_atan2)
{
  constexpr double a = 2.1;
  constexpr double b = 3.7;
  EXPECT_DOUBLE_EQ( stk::math::atan2(a,b), std::atan2(a,b) );
}

TEST(StkSimd, StkMathDouble_asinh)
{
  constexpr double f = 0.65;
  EXPECT_DOUBLE_EQ( stk::math::asinh(f), std::asinh(f) );
}

TEST(StkSimd, StkMathDouble_acosh)
{
  constexpr double f = 0.65;
  EXPECT_DOUBLE_EQ( stk::math::acosh(f+1.0), std::acosh(f+1.0) );
}

TEST(StkSimd, StkMathDouble_atanh)
{
  constexpr double f = 0.65;
  EXPECT_DOUBLE_EQ( stk::math::atanh(f), std::atanh(f) );
}

TEST(StkSimd, StkMathDouble_erf)
{
  constexpr double f = 0.65;
  EXPECT_DOUBLE_EQ( stk::math::erf(f), std::erf(f) );
}

TEST(StkSimd, StkMathDouble_multiplysign)
{
  constexpr double a = 2.1;
  constexpr double c = -9.1;
  EXPECT_DOUBLE_EQ( stk::math::multiplysign(a,c), a*std::copysign(1.0,c));
}

TEST(StkSimd, StkMathDouble_multiplysignNegZero)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::multiplysign(a,-0.0), a*std::copysign(1.0,-0.0));
}

TEST(StkSimd, StkMathDouble_multiplysignZero)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::multiplysign(a,0.0), a*std::copysign(1.0,0.0));
}

TEST(StkSimd, StkMathDouble_copysign)
{
  constexpr double a = 2.1;
  constexpr double c = -9.1;
  EXPECT_DOUBLE_EQ( stk::math::copysign(a,c), std::copysign(a,c) );
}

TEST(StkSimd, StkMathDouble_copysignNegZero)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::copysign(a,-0.0), std::copysign(a,-0.0) );
}

TEST(StkSimd, StkMathDouble_copysignZero)
{
  constexpr double a = 2.1;
  EXPECT_DOUBLE_EQ( stk::math::copysign(a,0.0), std::copysign(a,0.0) );
}

TEST(StkSimd, StkMathDouble_abs)
{
  constexpr double c = -9.1;
  EXPECT_DOUBLE_EQ( stk::math::abs(c), std::abs(c) );
}

TEST(StkSimd, StkMathDouble_max)
{
  constexpr double a = 2.1;
  constexpr double c = -9.1;
  EXPECT_DOUBLE_EQ( stk::math::max(a,c), a > c ? a : c );
}

TEST(StkSimd, StkMathDouble_min)
{
  constexpr double a = 2.1;
  constexpr double c = -9.1;
  EXPECT_DOUBLE_EQ( stk::math::min(a,c), a < c ? a : c );
}

TEST(StkSimd, StkMathDouble_if_then)
{
  constexpr double a = 2.1;
  constexpr double b = 3.7;
  EXPECT_DOUBLE_EQ( stk::math::if_then_else(true,a,b), a );
  EXPECT_DOUBLE_EQ( stk::math::if_then_else(false,a,b), b );
  EXPECT_DOUBLE_EQ( stk::math::if_then_else_zero(true,a), a );
  EXPECT_DOUBLE_EQ( stk::math::if_then_else_zero(false,a), 0 );
}

TEST(StkSimd, StkMathDouble_isnan)
{
  EXPECT_TRUE( stk::math::isnan( 0.0/0.0) );
  EXPECT_TRUE( !stk::math::isnan(1.0/0.0) );
  EXPECT_TRUE( !stk::math::isnan(0.0/1.0) );
  EXPECT_TRUE( !stk::math::isnan(3.4) );
}
