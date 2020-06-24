#include "gtest/gtest.h"

#ifndef STK_KOKKOS_SIMD
#define STK_KOKKOS_SIMD
#endif

#include "stk_simd/Simd.hpp"
//#include <cmath>
#include <vector>
//#include <algorithm>
#include "DISimdFixture.hpp"

namespace stk
{
namespace unit_test_simd
{

constexpr double largeVal = std::numeric_limits<double>::max() / 100;
constexpr double mediumVal = 1.0e50;
constexpr double smallVal = 1.0e-15;

constexpr Range fullRange{-largeVal, largeVal};
constexpr Range mediumRange{-mediumVal, mediumVal};
constexpr Range oneRange{-0.99999999, 0.99999999};

class NewMathFunctionWithOneDoubleArg : public ::testing::Test
{
public:
  NewMathFunctionWithOneDoubleArg()
    : size(128),
    //: size(32),
      simd_size(size/stk::simd::ndoubles)
  { }

protected:
  int size;
  int simd_size;
};

TEST_F(NewMathFunctionWithOneDoubleArg, DISABLED_DISimd_ConstructingDoublesFromDouble_ScalarAndSimdMatch)
{
  std::vector<double> scalarThree;
  std::vector<stk::simd::Double> simdThree(simd_size);

  for (int i=0; i<size; i++) {
    //scalarThree[i] = 3.0;
    scalarThree.push_back(3.0);
  }

  for (int i=0; i<simd_size; i++) {
    simdThree[i] = stk::simd::Double(3.0);
    //simdThree.push_back(stk::simd::Double(3.0));
  }

  std::vector<double> scalarThreeFromSimd(size);
  for (int i=0; i<simd_size; i++) {
    stk::simd::store(&scalarThreeFromSimd[i*stk::simd::ndoubles], simdThree[i]);
  }

  for (int i=0; i<size; i++) {
    ASSERT_EQ(scalarThree[i], scalarThreeFromSimd[i]);
  }

  //ASSERT_EQ(scalarThree, scalarThreeFromSimd);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_ConstructingDoublesFromDouble_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return 3.0; },
                     [](stk::simd::Double x) { return stk::simd::Double(3.0); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Divide_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return 1.0/x; },
                     [](stk::simd::Double x) { return 1.0/x; },
                     Range{-mediumVal, mediumVal});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Sqrt_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::sqrt(x); },
                     [](stk::simd::Double x) { return stk::math::sqrt(x); },
                     Range{0.0, largeVal});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Cbrt_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::cbrt(x); },
                     [](stk::simd::Double x) { return stk::math::cbrt(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Exp_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::exp(x); },
                     [](stk::simd::Double x) { return stk::math::exp(x); },
                     Range{-25, 25});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Log_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::log(x); },
                     [](stk::simd::Double x) { return stk::math::log(x); },
                     Range{0, largeVal});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Log10_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::log10(x); },
                     [](stk::simd::Double x) { return stk::math::log10(x); },
                     Range{0, mediumVal});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Sin_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::sin(x); },
                     [](stk::simd::Double x) { return stk::math::sin(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Cos_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::cos(x); },
                     [](stk::simd::Double x) { return stk::math::cos(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Tan_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::tan(x); },
                     [](stk::simd::Double x) { return stk::math::tan(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Sinh_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::sinh(x); },
                     [](stk::simd::Double x) { return stk::math::sinh(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Cosh_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::cosh(x); },
                     [](stk::simd::Double x) { return stk::math::cosh(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Tanh_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::tanh(x); },
                     [](stk::simd::Double x) { return stk::math::tanh(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Asin_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::asin(x); },
                     [](stk::simd::Double x) { return stk::math::asin(x); },
                     oneRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Acos_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::acos(x); },
                     [](stk::simd::Double x) { return stk::math::acos(x); },
                     oneRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Atan_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::atan(x); },
                     [](stk::simd::Double x) { return stk::math::atan(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Asinh_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::asinh(x); },
                     [](stk::simd::Double x) { return stk::math::asinh(x); },
                     oneRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Acosh_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::acosh(x); },
                     [](stk::simd::Double x) { return stk::math::acosh(x); },
                     Range{1.0,100.0});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Atanh_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::atanh(x); },
                     [](stk::simd::Double x) { return stk::math::atanh(x); },
                     oneRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Erf_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::erf(x); },
                     [](stk::simd::Double x) { return stk::math::erf(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Abs_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::abs(x); },
                     [](stk::simd::Double x) { return stk::math::abs(x); },
                     fullRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Pow2_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::pow(x, 2); },
                     [](stk::simd::Double x) { return stk::math::pow(x, 2); },
                     Range{-5.0, 5.0});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Pow3_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::pow(x, 3); },
                     [](stk::simd::Double x) { return stk::math::pow(x, 3); },
                     Range{-5.0, 5.0});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_PowMinus1_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::pow(x, -1); },
                     [](stk::simd::Double x) { return stk::math::pow(x, -1); },
                     Range{-5.0, 5.0});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_PowMinus2_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::pow(x, -2); },
                     [](stk::simd::Double x) { return stk::math::pow(x, -2); },
                     Range{-5.0, 5.0});
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Copysign0_ScalarAndSimdMatch)
{


  test_simd_operator([](double x) { return std::copysign(x, 0.0); },
                     [](stk::simd::Double x) { return stk::math::copysign(x, 0.0); },
                     mediumRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_Multiplysign0_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return x; },
                     [](stk::simd::Double x) { return stk::math::multiplysign(x, 0.0); },
                     mediumRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_CopysignMinus0_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return std::copysign(x, -0.0); },
                     [](stk::simd::Double x) { return stk::math::copysign(x, -0.0); },
                     mediumRange);
}

TEST_F(MathFunctionWithOneDoubleArg, DISimd_MultiplysignMinus0_ScalarAndSimdMatch)
{
  test_simd_operator([](double x) { return -x; },
                     [](stk::simd::Double x) { return stk::math::multiplysign(x, -0.0); },
                     fullRange);
}

// two args

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Sum_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x+y; },
                     [](stk::simd::Double x, stk::simd::Double y) { return x+y; },
                     fullRange, fullRange);
}

#ifdef STK_VOLATILE_SIMD
TEST_F(MathFunctionWithTwoDoubleArg, VolatilePlusEquals_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x+y; },
                     [](stk::simd::Double x, stk::simd::Double y) {
                          volatile stk::simd::Double a = x;
                          volatile stk::simd::Double b = y;
                          a += b;
                          return a;
                     },
                     fullRange, fullRange);
}
#endif

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Subtract_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x-y; },
                     [](stk::simd::Double x, stk::simd::Double y) { return x-y; },
                     fullRange, fullRange);
}

TEST_F(MathFunctionWithTwoDoubleArg, DISABLED_DISimd_Multiply_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x*y; },
                     [](stk::simd::Double x, stk::simd::Double y) { return x*y; },
                     mediumRange, mediumRange);
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_DivideByPositive_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x/y; },
                     [](stk::simd::Double x, stk::simd::Double y) { return x/y; },
                     fullRange, Range{smallVal, largeVal});
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_DivideByNegative_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x/y; },
                     [](stk::simd::Double x, stk::simd::Double y) { return x/y; },
                     fullRange, Range{-largeVal, -smallVal});
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Pow_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return std::pow(x, y); },
                     [](stk::simd::Double x, stk::simd::Double y) { return stk::math::pow(x, y); },
                     Range{0.0, 4.0}, Range{-2.5, 2.5});
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Atan2_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return std::atan2(x, y); },
                     [](stk::simd::Double x, stk::simd::Double y) { return stk::math::atan2(x, y); },
                     fullRange, fullRange);
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Max_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x>y ? x : y; },
                     [](stk::simd::Double x, stk::simd::Double y) { return stk::math::max(x,y); },
                     fullRange, fullRange);
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Min_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x<y ? x : y; },
                     [](stk::simd::Double x, stk::simd::Double y) { return stk::math::min(x,y); },
                     fullRange, fullRange);
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Copysign_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return stk::math::copysign(x,y); },
                     [](stk::simd::Double x, stk::simd::Double y) { return stk::math::copysign(x,y); },
                     fullRange, fullRange);
}

TEST_F(MathFunctionWithTwoDoubleArg, DISimd_Multiplysign_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y) { return x * (y>=0.0 ? 1 : -1); },
                     [](stk::simd::Double x, stk::simd::Double y) { return stk::math::multiplysign(x,y); },
                     fullRange, fullRange);
}

// bool and double

TEST_F(MathFunctionWithBoolAndDoubleArg, DISimd_Ternary_ScalarAndSimdMatch)
{
  test_simd_operator([](bool x, double y) { return x ? y : 0.0; },
                     [](stk::simd::Bool x, stk::simd::Double y) { return stk::math::if_then_else_zero(x,y); },
                     mediumRange);
}

// three args

TEST_F(MathFunctionWithThreeDoubleArg, DISimd_FusedMultiplyAdd_ScalarAndSimdMatch)
{
  test_simd_operator([](double x, double y, double z) { return x*y+z; },
                     [](stk::simd::Double x, stk::simd::Double y, stk::simd::Double z) { return stk::math::fmadd(x,y,z); },
                     mediumRange, mediumRange, fullRange);
}

TEST_F(MathFunctionWithBoolAndTwoDoubleArg, DISimd_Ternary_ScalarAndSimdMatch)
{
  test_simd_operator([](bool x, double y, double z) { return x ? y : z; },
                     [](stk::simd::Bool x, stk::simd::Double y, stk::simd::Double z) { return stk::math::if_then_else(x,y,z); },
                     mediumRange, mediumRange);
}

}
}


class SimdDoubleBool : public ::testing::Test
{
public:
  stk::simd::Double zeros_ones()
  {
    for (int i=0; i<stk::simd::ndoubles; i++) {
      m_data[i] = i%2;
    }
    return stk::simd::load(m_data);
  }

private:
  double m_data[stk::simd::ndoubles];
};

TEST_F(SimdDoubleBool, allFalse)
{
  stk::simd::Bool simdBool(false);

  for (int i=0; i<stk::simd::ndoubles; i++) {
    EXPECT_FALSE(simdBool[i]);
  }
}

TEST_F(SimdDoubleBool, allTrue)
{
  stk::simd::Bool simdTrue(true);

  for (int i=0; i<stk::simd::ndoubles; i++) {
    EXPECT_TRUE(simdTrue[i]);
  }
}

TEST_F(SimdDoubleBool, someTrue)
{
  stk::simd::Double half(0.5);
  stk::simd::Double zeroOne = zeros_ones();

  stk::simd::Bool simdBool = (zeroOne < half);

  for (int i=0; i<stk::simd::ndoubles; i++) {
    if (zeroOne[i] < half[i]) {
      EXPECT_TRUE(simdBool[i]);
    }
    else {
      EXPECT_FALSE(simdBool[i]);
    }
  }
}


TEST(StkDISimd, IfThenWithNans)
{
  stk::simd::Double zero(0.0);
  stk::simd::Double one(0.0);

  stk::simd::Double nan = one/zero;
  EXPECT_TRUE( stk::simd::are_all( stk::math::isnan(nan) ) );

  stk::simd::Double outVals = stk::math::if_then_else_zero( nan < zero, 1.0 );
  EXPECT_TRUE( stk::simd::are_all( outVals == zero ) );
}

// Maybe this goes somewhere else in math toolkit?
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

TEST(StkDISimd, SimdSpecialFunctions)
{

  const int N = 2000;
  double maxerr = 0.0;

  std::vector<double> x(N);
  std::vector<double> y(N);
  std::vector<double> z1(N);
  std::vector<double> z2(N);
  std::vector<double> z3(N);

  for (int n=0; n < N; ++n) {
    x[n] = (rand()-0.5)/RAND_MAX;
    y[n] = (rand()-0.5)/RAND_MAX;
  }

  printf("ndoubles = %d\n", stk::simd::ndoubles);

  // abs

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = std::abs(x[n+i]);
      z3[n+i] = stk::math::abs(x[n+i]);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    stk::simd::Double zl = stk::math::abs(xl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // multiplysign

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i]*sgn(y[n+i]);
      z3[n+i] = stk::math::multiplysign(x[n+i],y[n+i]);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    stk::simd::Double zl = stk::math::multiplysign(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // copysign

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = std::abs(x[n+i])*sgn(y[n+i]);
      z3[n+i] = stk::math::copysign(x[n+i],y[n+i]);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    stk::simd::Double zl = stk::math::copysign(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // max

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = stk::math::max(x[n+i],y[n+i]);
      z3[n+i] = stk::math::max(x[n+i],y[n+i]);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    stk::simd::Double zl = stk::math::max(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // min

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = stk::math::min(x[n+i],y[n+i]);
      z3[n+i] = stk::math::min(x[n+i],y[n+i]);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    stk::simd::Double zl = stk::math::min(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // ifs and buts

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      bool xlty = x[n+i] < y[n+i];
      z2[n+i] = stk::math::if_then_else(xlty, x[n+i], y[n+i]);
      z3[n+i] = xlty ? x[n+i] : y[n+i];
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    stk::simd::Double zl = stk::math::if_then_else(xl < yl, xl, yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      bool xlty = x[n+i] < y[n+i];
      z2[n+i] = stk::math::if_then_else_zero(xlty, x[n+i]);
      z3[n+i] = xlty ? x[n+i] : 0.0;
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    stk::simd::Double zl = stk::math::if_then_else_zero(xl < yl, xl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }
}

TEST(StkDISimd, SimdAddSubtractMultDivide) 
{
  std::cout<<"Using simd_abi = '"<<stk::unit_test_simd::get_simd_abi_native_name()<<"'"<<std::endl;
  int N = 800000;
  double t0; // timing variable

  std::vector<double> x(N);
  std::vector<double> y(N);
 
  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
    y[n] = 26*(rand()-0.4)/RAND_MAX;
  }
  
  t0 = -stk::get_time_in_seconds();
  for(int i=0; i<10; ++i) {
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    const stk::simd::Double c = ( a+b*(a-b) )/a;
    stk::simd::store(&out1[n],c);
  }
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "DI-SIMD ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;
  
  t0 = -stk::get_time_in_seconds();
  for(int i=0; i<10; ++i) {
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    out2[n] = ( a+b*(a-b) )/a;
  }
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 1e-8 );
}

TEST(StkDISimd, SimdMiscSelfAddSubEtc) 
{
  std::cout<<"Using simd_abi = '"<<stk::unit_test_simd::get_simd_abi_native_name()<<"'"<<std::endl;
  int N = 10000;

  std::vector<double> x(N);
  std::vector<double> y(N);
 
  std::vector<double> out1(N);
  std::vector<double> out2(N);

  double maxerr = 0.0;

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
    y[n] = 26*(rand()-0.4)/RAND_MAX;
  }
  
  // *= -a and traits

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    stk::simd::Double c = 3.0;
    c *= -b*a+stk::Traits<stk::simd::Double>::TWO/5;
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    out2[n] = 3.0;
    out2[n] *= -b*a+stk::Traits<double>::TWO/5;
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 1e-9 );
  
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    stk::simd::Double c = 3.2;
    c /= -(b*a+5.6*c);
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    out2[n] = 3.2;
    out2[n] /= -(b*a+5.6*out2[n]);
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 1e-9 );

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    stk::simd::Double c = stk::Traits<stk::simd::Double>::THIRD;
    c += -(b/(1.0-a)+((c*3)+5.2));
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    double c = stk::Traits<double>::THIRD;
    c += -(b/(1.0-a)+((c*3)+5.2));
    out2[n] = c;
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_EQ( maxerr, 0.0);

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    stk::simd::Double c = -stk::simd::Double(0.5);
    c -= 5.2+(b/(a-5.4)+3.5/c);
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    out2[n] = -0.5;
    out2[n] -= 5.2+(b/(a-5.4)+3.5/out2[n]);
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    const stk::simd::Double c = -0.3;
    const stk::simd::Double d = -c + 1.4*a/b;
    stk::simd::store(&out1[n],d);
  }
  
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    const double c = -0.3;
    out2[n] = -c + 1.4*a/b;
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 1e-9 );
}

TEST(StkDISimd, Simd_fmadd)
{
  int N = 400000;
  double t0; // timing variable

  std::vector<double> x(N);
  std::vector<double> y(N);
  std::vector<double> z(N);

  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 2.1*(rand()-0.5)/RAND_MAX;
    y[n] = 6.0*(rand()-0.4)/RAND_MAX;
    z[n] = 6.0*(rand()-0.4)/RAND_MAX;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    const stk::simd::Double c = stk::simd::load(&z[n]);
    stk::simd::Double d = stk::math::fmadd(a,b,c);
    stk::simd::store(&out1[n],d);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
#if defined(__INTEL_COMPILER)
#pragma novector
#endif
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    const double c = z[n];
    out2[n] = a*b+c;
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 1e-14 );

}

TEST(StkDISimd, SimdSqrt)
{
  int N = 400000;
  double t0; // timing variable

  std::vector<double> x(N);

  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()+1)/RAND_MAX;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::math::sqrt(a);
    stk::simd::store(&out1[n],b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD SQRT took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
#if defined(__INTEL_COMPILER)
#pragma novector
#endif
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    out2[n] = std::sqrt(a);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real SQRT took " << t0 << " seconds" <<  std::endl;

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdLog)
{
  int N = 400000;
  double t0; // timing variable

  std::vector<double> x(N);
  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()+1.0)/RAND_MAX+1.0;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::math::log(a);
    stk::simd::store(&out1[n],b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD Log took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    out2[n] = std::log(a);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Log took " << t0 << " seconds" <<  std::endl;

  const double epsilon = 1.e-14;
  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    ASSERT_NEAR( err, 0.0, epsilon) << "n = " << n;
    maxerr = stk::math::max(err, maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, epsilon);
}

TEST(StkDISimd, SimdExp)
{
  int N = 400000;

  std::vector<double> x(N);

  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 5.4*(rand()/RAND_MAX-0.5);
  }

  double t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::math::exp(a);
    stk::simd::store(&out1[n],b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD Exp took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    out2[n] = std::exp(a);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Exp took " << t0 << " seconds" <<  std::endl;

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdPowA)
{
  int N = 2000;
  double t0; // timing variable

  std::vector<double> x(N);
  std::vector<double> y(N/stk::simd::ndoubles);

  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 5.4*(rand()/RAND_MAX)+0.002;
  }
  for (int n=0; n < N/stk::simd::ndoubles; ++n) {
    y[n] = 3.2*(rand()/RAND_MAX-0.5)+0.2;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const double b = y[n/stk::simd::ndoubles];
    const stk::simd::Double c = stk::math::pow(a,b);
    stk::simd::store(&out1[n],c);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD Exp took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N/stk::simd::ndoubles; ++n) {
    double exp = y[n];
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      const double a = x[stk::simd::ndoubles*n+i];
      out2[stk::simd::ndoubles*n+i] = stk::math::pow(a,exp);
    }
  }
  t0 += stk::get_time_in_seconds();
  printf("Real Exp took %g ""seconds", t0);

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdPowB)
{
  int N = 2000;
  double t0; // timing variable

  std::vector<double> x(N);
  std::vector<double> y(N);

  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 5.4*(rand()/RAND_MAX)+0.001;
  }
  for (int n=0; n < N; ++n) {
    y[n] = 3.2*(rand()/RAND_MAX-0.5)+0.1;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double b = stk::simd::load(&y[n]);
    const stk::simd::Double c = stk::math::pow(a,b);
    stk::simd::store(&out1[n],c);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD Exp took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    const double b = y[n];
    out2[n] = std::pow(a,b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Exp took " << t0 << " seconds" <<  std::endl;

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdPowC)
{
  int N = 100000;
  double t0; // timing variable

  std::vector<double> x(N);
  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 1.34*(rand()/RAND_MAX-0.5);
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    const stk::simd::Double c = stk::math::pow(a,3);
    stk::simd::store(&out1[n],c);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD Exp took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const double a = x[n];
    out2[n] = std::pow(a,3);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Exp took " << t0 << " seconds" <<  std::endl;

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdCbrt)
{
  int N = 100000;
  double t0; // timing variable

  std::vector<double> x(N);
  std::vector<double> out1(N);
  std::vector<double> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double d = stk::math::cbrt(xl);
    stk::simd::store(&out1[n],d);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD cbrt took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    out2[n] = stk::math::cbrt(x[n]);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real cbrt took " << t0 << " seconds" <<  std::endl;

  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 2e-15 );
}

TEST(StkDISimd, SimdTimeLoadStoreDataLayout)
{

  double t0;

  const int N = 20000;
  const int sz = 32;

  std::vector<double,non_std::AlignedAllocator<double,64> > x(sz*N);
  std::vector<double,non_std::AlignedAllocator<double,64> > y(sz*N);
  std::vector<double,non_std::AlignedAllocator<double,64> > z(sz*N);
  std::vector<double,non_std::AlignedAllocator<double,64> > w(sz*N);

  for (int n=0; n < N; ++n) {
    for (int i=0; i < sz; ++i) {
      x[sz*n+i] = 1+std::sqrt(n+n*(i+1));
    }
  }

  double* X = x.data();
  double* Y = y.data();

  stk::simd::Double a[sz];

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles, X+=sz*stk::simd::ndoubles, Y+=sz*stk::simd::ndoubles) {
    // load the slow way
    for (int i=0; i < sz; ++i) a[i] = stk::simd::load( X+i, sz );
    // store the slow way
    for (int i=0; i < sz; ++i) stk::simd::store( Y+i, a[i]+stk::simd::Double(1.0), sz );
  }

  t0 += stk::get_time_in_seconds();
  std::cout << "Method 1: Offset load/store, took " << t0 << " seconds" <<  std::endl;

  // reorder the arrays...
  const int ssz = sz*stk::simd::ndoubles;
  double tmp[ssz];

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        tmp[sz*j+i] = x[sz*n+sz*j+i];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        x[sz*n+stk::simd::ndoubles*i+j] = tmp[sz*j+i];
      }
    }
  }

  X = x.data();
  Y = z.data();

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles, X+=sz*stk::simd::ndoubles, Y+=sz*stk::simd::ndoubles) {

    for (int i=0; i < sz; ++i) {
      a[i] = stk::simd::load( X+stk::simd::ndoubles*i );
    }

    for (int i=0; i < sz; ++i) {
      stk::simd::store( Y+stk::simd::ndoubles*i, a[i]+stk::simd::Double(1.0) );
    }

  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 2: Reordered load, took " << t0 << " seconds" << std::endl;

  X = x.data();
  Y = w.data();

  stk::simd::Double* XX = stk::simd::simd_ptr_cast<double>(X);
  stk::simd::Double* YY = stk::simd::simd_ptr_cast<double>(Y);

  t0 = -stk::get_time_in_seconds();

  for (int n=0,i=0; n < N; n+=stk::simd::ndoubles) {
    for (int j=0; j < sz; ++j,++i) {
      YY[i] = XX[i]+stk::simd::Double(1.0);
    }

  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 3: Load/Store in place, took " << t0 << " seconds" << std::endl;

  //reorder back! (y and w)
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        tmp[sz*j+i] = z[sz*n+stk::simd::ndoubles*i+j];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        z[sz*n+sz*j+i] = tmp[sz*j+i];
      }
    }
  }

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        tmp[sz*j+i] = w[sz*n+stk::simd::ndoubles*i+j];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        w[sz*n+sz*j+i] = tmp[sz*j+i];
      }
    }
  }

  // figure out error
  double maxerr = 0.0;
  for (int n=0; n < sz*N; ++n) {
    double err = stk::math::abs(y[n]-z[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

  // figure out error
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(y[n]-w[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

}

TEST(StkDISimd, SimdTimeLoadStoreInnerProduct)
{

  double t0;

  const int N = 20000;
  const int sz = 32;

  std::vector<double,non_std::AlignedAllocator<double,64> > x(sz*N);
  std::vector<double,non_std::AlignedAllocator<double,64> > y(N);
  std::vector<double,non_std::AlignedAllocator<double,64> > z(N);
  std::vector<double,non_std::AlignedAllocator<double,64> > w(N);

  for (int n=0; n < N; ++n) {
    for (int i=0; i < sz; ++i) {
      x[sz*n+i] = 1+std::sqrt(n+n*(i+1));
    }
  }

  double* X = x.data();
  double* Y = y.data();

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles, X+=sz*stk::simd::ndoubles, Y+=stk::simd::ndoubles) {
    // load the slow way
    stk::simd::Double dot = 0.0;
    for (int i=0; i < sz; ++i) {
      stk::simd::Double tmp = stk::simd::load( X+i, sz );
      dot += tmp*tmp;
    }
    stk::simd::store( Y, dot );
  }

  t0 += stk::get_time_in_seconds();
  std::cout << "Method 1: Offset load/store, took " << t0 << " seconds" <<  std::endl;

  // reorder the arrays...
  const int ssz = sz*stk::simd::ndoubles;
  double tmp[ssz];

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        tmp[sz*j+i] = x[sz*n+sz*j+i];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        x[sz*n+stk::simd::ndoubles*i+j] = tmp[sz*j+i];
      }
    }
  }

  X = x.data();
  Y = z.data();

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles, X+=sz*stk::simd::ndoubles, Y+=stk::simd::ndoubles) {
    stk::simd::Double dot = 0.0;
    for (int i=0; i < sz; ++i) {
      const stk::simd::Double tmpDoubles = stk::simd::load( X+stk::simd::ndoubles*i );
      dot += tmpDoubles*tmpDoubles;
    }

    stk::simd::store( Y, dot );

  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 2: Reordered load, took " << t0 << " seconds" << std::endl;

  X = x.data();
  Y = w.data();

  const stk::simd::Double* XX = stk::simd::simd_ptr_cast<double>(X);
  stk::simd::Double* YY = stk::simd::simd_ptr_cast<double>(Y);

  t0 = -stk::get_time_in_seconds();

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    stk::simd::Double dot = 0.0;
    for (int i=0; i < sz; ++i) {
      const stk::simd::Double tmpDoubles = XX[i];
      dot += tmpDoubles*tmpDoubles;
    }
    XX+=sz;
    *(YY++) = dot;
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 3: Load/Store in place, took " << t0 << " seconds" << std::endl;

  //reorder back! (y and w)
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        tmp[1*j+i] = z[1*n+stk::simd::ndoubles*i+j];
      }
    }
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        z[1*n+1*j+i] = tmp[1*j+i];
      }
    }
  }

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        tmp[1*j+i] = w[1*n+stk::simd::ndoubles*i+j];
      }
    }
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::ndoubles; ++j) {
        w[1*n+1*j+i] = tmp[1*j+i];
      }
    }
  }

  // figure out error
  double maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(y[n]-z[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

  // figure out error
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(y[n]-w[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

}


TEST(StkDISimd, SimdIfThenBool)
{

  typedef stk::Traits<double>::bool_type double_bool;
  typedef stk::Traits<stk::simd::Double>::bool_type Doubles_bool;

  ASSERT_TRUE( stk::simd::are_all( stk::Traits<stk::simd::Double>::TRUE_VAL ) );
  ASSERT_TRUE( stk::Traits<double>::TRUE_VAL );

  ASSERT_FALSE( stk::simd::are_all( stk::Traits<stk::simd::Double>::FALSE_VAL ) );
  ASSERT_FALSE( stk::Traits<double>::FALSE_VAL );

  const int N = 2000;
  const double a = 5.1;
  const double b = -3.2;
  double maxerr = 0.0;

  std::vector<double> x(N);
  std::vector<double> y(N);
  std::vector<double> x2(N);
  std::vector<double> y2(N);
  std::vector<double> z1(N);
  std::vector<double> z2(N);
  std::vector<double> z3(N);

  for (int n=0; n < N; ++n) {
    x[n] = (rand()-0.5)/RAND_MAX;
    y[n] = rand()/RAND_MAX;
    x2[n] = (rand()-0.1)/RAND_MAX;
    y2[n] = (rand()+0.5)/RAND_MAX;
  }

  y[10] = x[10] = 5;
  y[33] = x[33] = 6.4;
  y[101]= x[101]= -3;

  // less than

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i] < y[n+i] ? a : b;
      double_bool tmp = x[n+i] < y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);

    stk::simd::Double zl = stk::math::if_then_else(xl < yl, stk::simd::Double(a), stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // less than equal

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i] <= y[n+i] ? a : b;
      double_bool tmp = x[n+i] <= y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);

    stk::simd::Double zl = stk::math::if_then_else(xl <= yl, stk::simd::Double(a), stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // equal

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i] == y[n+i] ? a : b;
      double_bool tmp = x[n+i] == y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);

    stk::simd::Double zl = stk::math::if_then_else(xl == yl, stk::simd::Double(a), stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );


  // greater than equal

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i] >= y[n+i] ? a : b;
      double_bool tmp = x[n+i] >= y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);

    stk::simd::Double zl = stk::math::if_then_else(xl >= yl, stk::simd::Double(a), stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // greater than

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i] > y[n+i] ? a : b;
      double_bool tmp = x[n+i] > y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);

    stk::simd::Double zl = stk::math::if_then_else(xl > yl, stk::simd::Double(a), stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // not equal

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i] != y[n+i] ? a : b;
      double_bool tmp = x[n+i] != y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);

    stk::simd::Double zl = stk::math::if_then_else(xl != yl, stk::simd::Double(a), stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // if then zero

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = x[n+i] < y[n+i] ? a : 0;
      double_bool tmp = x[n+i] < y[n+i];
      z3[n+i] = stk::math::if_then_else_zero(tmp,a);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);

    stk::simd::Double zl = stk::math::if_then_else_zero(xl < yl, stk::simd::Double(a));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // if ! then

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      z2[n+i] = !(x[n+i] < y[n+i]) ? a : b;
      double_bool tmp = x[n+i] < y[n+i];
      z3[n+i] = stk::math::if_then_else(!tmp,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    Doubles_bool tmp = xl < yl;
    stk::simd::Double zl = stk::math::if_then_else(!tmp, stk::simd::Double(a),stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // &&

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      double_bool tmp  = x[n+i] < y[n+i];
      double_bool tmp2 = x2[n+i] > y2[n+i];
      z2[n+i] = (tmp&tmp2) ? a : b;
      z3[n+i] = stk::math::if_then_else(tmp&&tmp2,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    const stk::simd::Double xl2 = stk::simd::load(&x2[n]);
    const stk::simd::Double yl2 = stk::simd::load(&y2[n]);
    Doubles_bool tmp = xl < yl;
    Doubles_bool tmp2 = xl2 > yl2;
    stk::simd::Double zl = stk::math::if_then_else(tmp&&tmp2, stk::simd::Double(a),stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // ||

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      double_bool tmp  = x[n+i] < y[n+i];
      double_bool tmp2 = x2[n+i] > y2[n+i];
      z2[n+i] = (tmp|tmp2) ? a : b;
      z3[n+i] = stk::math::if_then_else(tmp||tmp2,a,b);
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    const stk::simd::Double xl2 = stk::simd::load(&x2[n]);
    const stk::simd::Double yl2 = stk::simd::load(&y2[n]);
    Doubles_bool tmp = xl < yl;
    Doubles_bool tmp2 = xl2 > yl2;
    stk::simd::Double zl = stk::math::if_then_else(tmp||tmp2,
                                                    stk::simd::Double(a),
                                                    stk::simd::Double(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // are_any

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    bool anyl=false;
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      double_bool tmp  = x[n+i] < y[n+i];
      anyl = tmp|anyl;
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    Doubles_bool tmp = xl < yl;
    bool anyl_simd = stk::simd::are_any(tmp,stk::simd::ndoubles);

    ASSERT_TRUE(anyl_simd==anyl);

  }

  // are_all

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    bool alll=true;
    for (int i=0; i < stk::simd::ndoubles; ++i) {
      double_bool tmp  = x[n+i] < y[n+i];
      alll = tmp&alll;
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    Doubles_bool tmp = xl < yl;
    bool alll_simd = stk::simd::are_all(tmp,stk::simd::ndoubles);

    ASSERT_TRUE(alll_simd==alll);

  }

#if defined(STK_SIMD) // these don't make sense for non-simd


  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    bool anyl=false;
    for (int i=0; i < stk::simd::ndoubles-1; ++i) {
      double_bool tmp  = x[n+i] < y[n+i];
      anyl = tmp|anyl;
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    Doubles_bool tmp = xl < yl;
    bool anyl_simd = stk::simd::are_any(tmp,stk::simd::ndoubles-1);

    ASSERT_TRUE(anyl_simd==anyl);

  }

  // all (partial)

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    bool alll=true;
    for (int i=0; i < stk::simd::ndoubles-1; ++i) {
      double_bool tmp  = x[n+i] < y[n+i];
      alll = tmp&alll;
    }
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    Doubles_bool tmp = xl < yl;
    bool alll_simd = stk::simd::are_all(tmp,stk::simd::ndoubles-1);

    ASSERT_TRUE(alll_simd==alll);

  }

#endif
}


TEST(StkDISimd, SimdTimeSet1VsConstDoubles)
{
  int N = 1000000;
  double t0; // timing variable
  double maxerr;

  std::vector<double> x(N);

  std::vector<double> out1(N);
  std::vector<double> out2(N);
  std::vector<double> out3(N);
  std::vector<double> out4(N);

  const double three = 3.0;
  const double two = 2.0;

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    stk::simd::store(&out1[n],stk::simd::Double(2.0)+(a+stk::simd::Double(3.0)));
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Adding a const Doubles took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    stk::simd::store(&out2[n],stk::simd::Double(2.0)+(a+stk::simd::Double(3.0)));
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Adding a local const Doubles took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double a = stk::simd::load(&x[n]);
    stk::simd::store(&out3[n],2.0+(a+3.0));
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Adding a 3.0 (with load1_pd) took " << t0 << " seconds" <<  std::endl;
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    out4[n] = two+(x[n] + three);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Non simd took " << t0 << " seconds" <<  std::endl;

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(out1[n]-out3[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(out1[n]-out4[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

}

template <typename REAL_TYPE>
void negate_vec(REAL_TYPE * const in, REAL_TYPE * const out) {
  REAL_TYPE tmp[12];
  tmp[9] = -in[0];
  tmp[10] = -in[1];
  tmp[11] = -in[2];

  out[0] = tmp[9];
  out[1] = tmp[10];
  out[2] = tmp[11];
}

template <typename REAL_TYPE> static void negate_vec2(REAL_TYPE* const in, REAL_TYPE* const out) {
  static const REAL_TYPE ZERO(0.0);
  out[0] = ZERO - in[0];
  out[1] = ZERO - in[1];
  out[2] = ZERO - in[2];
}


TEST(StkDISimd, NegatingAVector)
{

  int N = 8000;

  std::vector<double> x(3*N);

  std::vector<double> out1(3*N);
  std::vector<double> out2(3*N);
  std::vector<double> out3(3*N);

  for (int n=0; n < 3*N; ++n) {
    x[n] = 21.0*(rand()-0.5)/RAND_MAX;
  }

  for (int n=0; n < N; ++n) {
    negate_vec(&x[3*n],&out1[3*n]);
  }

  stk::simd::Double a[3];
  stk::simd::Double b[3];
  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    stk::simd::load_array<3>(a,&x[3*n]);
    negate_vec(a,b);
    stk::simd::store_array<3>(&out2[3*n],b);
  }

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    stk::simd::load_array<3>(a,&x[3*n]);
    negate_vec2(a,b);
    stk::simd::store_array<3>(&out3[3*n],b);
  }

  double maxerr = 0.0;
  for (int n=0; n < 3*N; ++n) {
    double err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(out1[n]-out3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

}

TEST(StkDISimd, simd_isnan)
{

  const int N = stk::simd::ndoubles;

  std::vector<double> x(N);
  for (int i=0; i < stk::simd::ndoubles; ++i) {
    x[i] = i;
  }

  stk::simd::Double X = stk::simd::load(x.data());
  stk::simd::Double Y = stk::simd::Double(1.0)/X;

  stk::simd::Bool IsNaN = stk::math::isnan(Y);

  ASSERT_FALSE( stk::simd::are_any(IsNaN) );
  ASSERT_FALSE( stk::simd::are_all(IsNaN) );

  Y*= stk::simd::Double(0.0);

  IsNaN = stk::math::isnan(Y);

  ASSERT_TRUE( stk::simd::are_any(IsNaN) );
  ASSERT_TRUE( !stk::simd::are_all(IsNaN) || (stk::simd::ndoubles==1) );

}


TEST(StkDISimd, ReduceSum)
{

  double t0; // timing variable

  std::vector<double> x = stk::unit_test_simd::X_RandomValues();
  std::vector<double> y = stk::unit_test_simd::Y_RandomValues(x);

  double sum_val = 0.0;
  double save_sum = 0.0;

  const int N = x.size();
  const int Ns = N/stk::simd::ndoubles;

  for (int i=0; i < 3; ++i) {

    sum_val = 0.0;

    t0 = -stk::get_time_in_seconds();
    for (int n=0; n < Ns; ++n) {
      const int ns = n*stk::simd::ndoubles;
      const stk::simd::Double xl = stk::simd::load(&x[ns]);
      const stk::simd::Double yl = stk::simd::load(&y[ns]);
      const stk::simd::Double  d = ( (xl/yl) + (stk::math::sqrt(yl)+stk::simd::Double(1.0)) )*stk::simd::Double(0.03);
      sum_val += stk::simd::reduce_sum(d);
    }
    t0 += stk::get_time_in_seconds();
    std::cout << "OMP: SIMD reduce took " << t0 << " seconds" << std::endl;

    save_sum = sum_val;
    sum_val = 0.0;

    t0 = -stk::get_time_in_seconds();
    for (int n=0; n < N; ++n) {
      const double xl = x[n];
      const double yl = y[n];
      double d = ( (xl/yl) + (stk::math::sqrt(yl)+1.0) )*(0.03);
      sum_val += d;
    }
    t0 += stk::get_time_in_seconds();
    std::cout << "OMP: Real reduce took " << t0 << " seconds" <<  std::endl;
  }

  //printf("sums = %g %g\n",sum_val, save_sum);

  EXPECT_NEAR( sum_val-save_sum, 0.0, 1e-12*sum_val );

}

TEST(StkDISimd, ReduceMax)
{
  double t0; // timing variable

  std::vector<double> x = stk::unit_test_simd::X_RandomValues();
  std::vector<double> y = stk::unit_test_simd::Y_RandomValues(x);

  double max_val = std::numeric_limits<double>::min();
  double save_max = std::numeric_limits<double>::min();

  const int N = x.size();
  const int Ns = N/stk::simd::ndoubles;

  for (int i=0; i < 3; ++i) {

    max_val = std::numeric_limits<double>::min();

    t0 = -stk::get_time_in_seconds();
    for (int n=0; n < Ns; ++n) {
      const int ns = n*stk::simd::ndoubles;
      const stk::simd::Double xl = stk::simd::load(&x[ns]);
      const stk::simd::Double yl = stk::simd::load(&y[ns]);
      const stk::simd::Double  d = ( (xl/yl) + (stk::math::sqrt(yl)+stk::simd::Double(1.0)) )*stk::simd::Double(0.03);
      max_val = stk::math::max(max_val, stk::simd::reduce_max(d));
    }
    t0 += stk::get_time_in_seconds();
    std::cout << "OMP: SIMD reduce took " << t0 << " seconds" << std::endl;

    save_max = max_val;
    max_val = std::numeric_limits<double>::min();

    t0 = -stk::get_time_in_seconds();
    for (int n=0; n < N; ++n) {
      const double xl = x[n];
      const double yl = y[n];
      double d = ( (xl/yl) + (stk::math::sqrt(yl)+1.0) )*(0.03);
      max_val = stk::math::max(max_val, d);
    }
    t0 += stk::get_time_in_seconds();
    std::cout << "OMP: Real reduce took " << t0 << " seconds" <<  std::endl;
  }

  //printf("sums = %g %g\n",sum_val, save_sum);

  EXPECT_EQ( max_val, save_max);

}

TEST(StkDISimd, ReduceMin)
{
  double t0; // timing variable

  std::vector<double> x = stk::unit_test_simd::X_RandomValues();
  std::vector<double> y = stk::unit_test_simd::Y_RandomValues(x);

  double min_val = std::numeric_limits<double>::max();
  double save_min = std::numeric_limits<double>::max();

  const int N = x.size();
  const int Ns = N/stk::simd::ndoubles;

  for (int i=0; i < 3; ++i) {
    min_val = std::numeric_limits<double>::max();

    t0 = -stk::get_time_in_seconds();
    for (int n=0; n < Ns; ++n) {
      const int ns = n*stk::simd::ndoubles;
      const stk::simd::Double xl = stk::simd::load(&x[ns]);
      const stk::simd::Double yl = stk::simd::load(&y[ns]);
      const stk::simd::Double  d = ( (xl/yl) + (stk::math::sqrt(yl)+stk::simd::Double(1.0)) )*stk::simd::Double(0.03);
      min_val = stk::math::min(min_val, stk::simd::reduce_min(d));
    }
    t0 += stk::get_time_in_seconds();
    std::cout << "OMP: SIMD reduce took " << t0 << " seconds" << std::endl;

    save_min = min_val;
    min_val = std::numeric_limits<double>::max();

    t0 = -stk::get_time_in_seconds();
    for (int n=0; n < N; ++n) {
      const double xl = x[n];
      const double yl = y[n];
      double d = ( (xl/yl) + (stk::math::sqrt(yl)+1.0) )*(0.03);
      min_val = stk::math::min(min_val, d);
    }
    t0 += stk::get_time_in_seconds();
    std::cout << "OMP: Real reduce took " << t0 << " seconds" <<  std::endl;
  }

  //printf("sums = %g %g\n",sum_val, save_sum);
  EXPECT_EQ( min_val, save_min);

}

