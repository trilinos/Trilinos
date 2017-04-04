#include <gtest/gtest.h>
#include <stk_math/StkMath.hpp>
#include <cmath>
#include <algorithm>

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
  EXPECT_NEAR( stk::math::cbrt(a), std::pow(a, float(1.0/3.0) ), 4e-7 );
}

TEST(StkSimd, StkMathFloat_log)
{
  EXPECT_EQ( stk::math::log(a), std::log(a) );
}

TEST(StkSimd, StkMathFloat_exp)
{
  EXPECT_EQ( stk::math::exp(a), std::exp(a) );
}

TEST(StkSimd, StkMathFloat_pow)
{
  EXPECT_NEAR( stk::math::pow(a,d), std::pow(a,d), 4e-7 );
  EXPECT_NEAR( stk::math::pow(a,b), std::pow(a,b), 4e-7 );
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

TEST(StkSimd, StkMathFloat_multiplysign)
{
  EXPECT_EQ( stk::math::multiplysign(a,c), a * (c>=0 ? 1.0 : -1.0) );
}

TEST(StkSimd, StkMathFloat_copysign)
{
  EXPECT_EQ( stk::math::copysign(a,c), std::copysign(a,c) );
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

TEST(StkSimd, StkMathFloat_if_not_then)
{
  EXPECT_EQ( stk::math::if_not_then_else(false,a,b), a );
  EXPECT_EQ( stk::math::if_not_then_else(true,a,b), b );
  EXPECT_EQ( stk::math::if_not_then_else_zero(false,a), a );
  EXPECT_EQ( stk::math::if_not_then_else_zero(true,a), 0.0f );
}

TEST(StkSimd, StkMathFloat_isnan)
{
  EXPECT_TRUE( stk::math::isnan(0.0f/0.0f) );
  EXPECT_TRUE( !stk::math::isnan(1.0f/0.0f) );
  EXPECT_TRUE( !stk::math::isnan(0.0f/1.0f) );
  EXPECT_TRUE( !stk::math::isnan(3.4f) );
}
