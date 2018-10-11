#include <gtest/gtest.h>
#include <stk_math/StkMath.hpp>
#include <cmath>
#include <algorithm>
#include <limits>

constexpr float a = 2.1f;
constexpr float b = 3.7f;
constexpr float c = -9.1f;
constexpr float f = 0.65f;
constexpr int d = 2;

TEST(StkSimd, StkMathFloat_fmadd)
{
  EXPECT_EQ( stk::math::fmadd(a,b,c), (a*b)+c );
}

TEST(StkSimd, StkMathFloat_sqrt)
{
  EXPECT_EQ( stk::math::sqrt(a), std::sqrt(a) );
}

TEST(StkSimd, StkMathFloat_cbrt)
{
  const float epsilon = 2 * std::numeric_limits<float>::epsilon();
  EXPECT_NEAR( stk::math::cbrt(a), std::pow(a, float(1.0/3.0) ), epsilon );
}

TEST(StkSimd, StkMathFloat_log)
{
  EXPECT_EQ( stk::math::log(a), std::log(a) );
}

TEST(StkSimd, StkMathFloat_log10)
{
  EXPECT_EQ( stk::math::log10(a), std::log10(a) );
}

TEST(StkSimd, StkMathFloat_exp)
{
  EXPECT_EQ( stk::math::exp(a), std::exp(a) );
}

TEST(StkSimd, StkMathFloat_pow)
{
  const float epsilon = 2 * std::numeric_limits<float>::epsilon();
  EXPECT_NEAR( stk::math::pow(a,d), std::pow(a,d), epsilon );
  EXPECT_NEAR( stk::math::pow(a,b), std::pow(a,b), epsilon );
}

TEST(StkSimd, StkMathFloat_sin)
{
  EXPECT_EQ( stk::math::sin(a), std::sin(a) );
}

TEST(StkSimd, StkMathFloat_cos)
{
  EXPECT_EQ( stk::math::cos(a), std::cos(a) );
}

TEST(StkSimd, StkMathFloat_tan)
{
  EXPECT_EQ( stk::math::tan(a), std::tan(a) );
}

TEST(StkSimd, StkMathFloat_sinh)
{
  EXPECT_EQ( stk::math::sinh(a), std::sinh(a) );
}

TEST(StkSimd, StkMathFloat_cosh)
{
  EXPECT_EQ( stk::math::cosh(a), std::cosh(a) );
}

TEST(StkSimd, StkMathFloat_tanh)
{
  EXPECT_EQ( stk::math::tanh(a), std::tanh(a) );
}

TEST(StkSimd, StkMathFloat_asin)
{
  EXPECT_EQ( stk::math::asin(f), std::asin(f) );
}

TEST(StkSimd, StkMathFloat_acos)
{
  EXPECT_EQ( stk::math::acos(f), std::acos(f) );
}

TEST(StkSimd, StkMathFloat_atan)
{
  EXPECT_EQ( stk::math::atan(f), std::atan(f) );
}

TEST(StkSimd, StkMathFloat_atan2)
{
  EXPECT_EQ( stk::math::atan2(a,b), std::atan2(a,b) );
}

TEST(StkSimd, StkMathFloat_asinh)
{
  EXPECT_EQ( stk::math::asinh(f), std::asinh(f) );
}

TEST(StkSimd, StkMathFloat_acosh)
{
  EXPECT_EQ( stk::math::acosh(f+1.f), std::acosh(f+1.f) );
}

TEST(StkSimd, StkMathFloat_atanh)
{
  EXPECT_EQ( stk::math::atanh(f), std::atanh(f) );
}

TEST(StkSimd, StkMathFloat_erf)
{
  EXPECT_EQ( stk::math::erf(f), std::erf(f) );
}

TEST(StkSimd, StkMathFloat_multiplysign)
{
  EXPECT_EQ( stk::math::multiplysign(a,c), a*std::copysign(1.0,c));
}

TEST(StkSimd, StkMathFloat_multiplysignNegZero)
{
  EXPECT_EQ( stk::math::multiplysign(a,-0.f), a*std::copysign(1.0,-0.f) );
}

TEST(StkSimd, StkMathFloat_multiplysignZero)
{
  EXPECT_EQ( stk::math::multiplysign(a,0.f), a*std::copysign(1.0,0.f) );
}

TEST(StkSimd, StkMathFloat_copysign)
{
  EXPECT_EQ( stk::math::copysign(a,c), std::copysign(a,c) );
}

TEST(StkSimd, StkMathFloat_copysignNegZero)
{
  EXPECT_EQ( stk::math::copysign(a,-0.f), std::copysign(a,-0.f) );
}

TEST(StkSimd, StkMathFloat_copysignZero)
{
  EXPECT_EQ( stk::math::copysign(a,0.f), std::copysign(a,0.f) );
}

TEST(StkSimd, StkMathFloat_abs)
{
  EXPECT_EQ( stk::math::abs(c), std::abs(c) );
}

TEST(StkSimd, StkMathFloat_max)
{
  EXPECT_EQ( stk::math::max(a,c), a > c ? a : c );
}

TEST(StkSimd, StkMathFloat_min)
{
  EXPECT_EQ( stk::math::min(a,c), a < c ? a : c );
}

TEST(StkSimd, StkMathFloat_if_then)
{
  EXPECT_EQ( stk::math::if_then_else(true,a,b), a );
  EXPECT_EQ( stk::math::if_then_else(false,a,b), b );
  EXPECT_EQ( stk::math::if_then_else_zero(true,a), a );
  EXPECT_EQ( stk::math::if_then_else_zero(false,a), 0.0f );
}

TEST(StkSimd, StkMathFloat_isnan)
{
  EXPECT_TRUE( stk::math::isnan(0.0f/0.0f) );
  EXPECT_TRUE( !stk::math::isnan(1.0f/0.0f) );
  EXPECT_TRUE( !stk::math::isnan(0.0f/1.0f) );
  EXPECT_TRUE( !stk::math::isnan(3.4f) );
}
