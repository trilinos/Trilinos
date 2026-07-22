// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef TAYLORUNITTESTS_HPP
#define TAYLORUNITTESTS_HPP

// ADOL-C includes
#include "adolc/adouble.h"
#include "adolc/interfaces.h"
#include "adolc/taping.h"

// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Random.hpp"

inline adouble max(const adouble& a, const adouble& b) { return fmax(a,b); }
inline adouble max(const adouble& a, double v) { return fmax(a,v); }
inline adouble max(double v, const adouble& b) { return fmax(v,b); }
inline adouble min(const adouble& a, const adouble& b) { return fmin(a,b); }
inline adouble min(const adouble& a, double v) { return fmin(a,v); }
inline adouble min(double v, const adouble& b) { return fmin(v,b); }

// gtest includes
#include <gtest/gtest.h>

#include "GTestUtils.hpp"

#define COMPARE_POLYS(x_dtay, x_adolc)                          \
  ASSERT_TRUE(x_dtay.degree() == d);                            \
  for (int i=0; i<=d; i++) {                                    \
    COMPARE_VALUES(x_dtay.coeff(i), x_adolc[i]);               \
  }                                                             \
  ;

#define COMPARE_TAYS(x_dtay, y_dtay)                            \
  ASSERT_TRUE(x_dtay.degree() == y_dtay.degree());              \
  for (int i=0; i<=x_dtay.degree(); i++) {                      \
    COMPARE_VALUES(x_dtay.coeff(i), y_dtay.coeff(i));          \
  }                                                             \
  ;

// A class for testing each Taylor operation
template <class TaylorType>
class TaylorOpsUnitTest : public ::testing::Test {
protected:

  // Taylor variables
  TaylorType a_dtay_, b_dtay_, c_dtay_;

  // ADOL-C arrays
  double **X_, **Y_;

  // Random number generator
  Sacado::Random<double> urand;

  // Degree of polynomials
  int d_;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

  TaylorOpsUnitTest()  :
    urand(), d_(5), tol_a(1.0e-11), tol_r(1.0e-10)
  {
    X_ = new double*[2];
    X_[0] = new double[d_+1];
    X_[1] = new double[d_+1];

    Y_ = new double*[1];
    Y_[0] = new double[d_+1];
  }

  ~TaylorOpsUnitTest()
  {
    delete [] X_[1];
    delete [] X_[0];
    delete [] X_;

    delete [] Y_[0];
    delete [] Y_;
  }

  void SetUp() {
    double val;

    a_dtay_ = TaylorType(d_,0.0);
    b_dtay_ = TaylorType(d_,0.0);

    for (int i=0; i<=d_; i++) {
      val = urand.number();
      a_dtay_.fastAccessCoeff(i) = val;
      X_[0][i] = val;

      val = urand.number();
      b_dtay_.fastAccessCoeff(i) = val;
      X_[1][i] = val;

      Y_[0][i] = 0.0;
    }
  }

  void TearDown() {}

  template <typename ScalarT>
  ScalarT composite1(const ScalarT& a, const ScalarT& b) {
    ScalarT t1 = 3. * a + sin(b) / log(fabs(a - b * 7.));
    ScalarT t2 = 1.0e3;
    ScalarT t3 = 5.7e4;
    ScalarT t4 = 3.2e5;
    t1 *= cos(a + exp(t1)) / 6. - tan(t1*sqrt(fabs(a * log10(fabs(b)))));
    t1 -= acos((6.+asin(pow(fabs(a),b)/t2))/t3) * asin(pow(fabs(b),2.)*1.0/t4) * atan((b*pow(2.,log(fabs(a))))/(t3*t4));
    t1 /= cosh(b - 0.7) + 7.*sinh(t1 + 0.8)*tanh(9./(a+1.)) - 9.;
    t1 += pow(fabs(a*4.),b-8.)/cos(a*b*a);

    return t1;
  }

  // void print_poly(double *x) {
  //   std::cout.setf(std::ios::fixed,std::ios::floatfield);
  //   std::cout.width(12);
  //   std::cout << "[";

  //   for (int i=0; i<=d_; i++) {
  //     std::cout.width(12);
  //     std::cout << x[i];
  //   }

  //   std::cout << "]\n";
  // }

  // void print_diff(const TaylorType& x_dtay, double* x_adolc)  {
  //   std::cout.setf(std::ios::fixed,std::ios::floatfield);
  //   std::cout.width(12);
  //   std::cout << "[";

  //   for (int i=0; i<=d_; i++) {
  //     std::cout.width(12);
  //     std::cout << x[i];
  //   }

  //   std::cout << "]\n";
  // }

}; // class TaylorOpsUnitTest

// Additional test fixture of max/min since CacheTaylor doesn't implement them
template <typename TaylorType>
class TaylorMaxMinUnitTest : public TaylorOpsUnitTest< TaylorType > {
protected:
  TaylorMaxMinUnitTest() {}
  ~TaylorMaxMinUnitTest() {}
};

TYPED_TEST_SUITE_P(TaylorOpsUnitTest);
TYPED_TEST_SUITE_P(TaylorMaxMinUnitTest);

#define BINARY_OP2_TEST(FIXTURENAME,TESTNAME,OP)        \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {  \
    auto a_dtay = this->a_dtay_;            \
    auto b_dtay = this->b_dtay_;            \
    auto c_dtay = this->c_dtay_;            \
    auto d = this->d_;                      \
    auto X = this->X_;                      \
    auto Y = this->Y_;                      \
    c_dtay = a_dtay OP b_dtay;              \
    trace_on(0);                            \
    adouble aa, ab, ac;                     \
    aa <<= X[0][0];                         \
    ab <<= X[1][0];                         \
    ac = aa OP ab;                          \
    ac >>= Y[0][0];                         \
    trace_off();                            \
    forward(0,1,2,d,0,X,Y);                 \
    COMPARE_POLYS(c_dtay,Y[0]);             \
  }

#define BINARY_OPRC_TEST(FIXTURENAME,TESTNAME,OP)       \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {   \
    auto a_dtay = this->a_dtay_;            \
    auto c_dtay = this->c_dtay_;            \
    auto d = this->d_;                      \
    auto X = this->X_;                      \
    auto Y = this->Y_;                      \
    double val = this->urand.number();      \
    c_dtay = a_dtay OP val;                 \
    trace_on(0);                            \
    adouble aa, ac;                         \
    aa <<= X[0][0];                         \
    ac = aa OP val;                         \
    ac >>= Y[0][0];                         \
    trace_off();                            \
    forward(0,1,1,d,0,X,Y);                 \
    COMPARE_POLYS(c_dtay,Y[0]);             \
  }

#define BINARY_OPLC_TEST(FIXTURENAME,TESTNAME,OP)       \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {   \
    auto a_dtay = this->a_dtay_;            \
    auto c_dtay = this->c_dtay_;            \
    auto d = this->d_;                      \
    auto X = this->X_;                      \
    auto Y = this->Y_;                      \
    double val = this->urand.number();      \
    c_dtay = val OP a_dtay;                 \
    trace_on(0);                            \
    adouble aa, ac;                         \
    aa <<= X[0][0];                         \
    ac = val OP aa;                         \
    ac >>= Y[0][0];                         \
    trace_off();                            \
    forward(0,1,1,d,0,X,Y);                 \
    COMPARE_POLYS(c_dtay,Y[0]);             \
  }

#define BINARY_OP_TEST(FIXTURENAME,TESTNAME,OP)                     \
  BINARY_OP2_TEST(FIXTURENAME,TESTNAME,OP)                          \
  BINARY_OPLC_TEST(FIXTURENAME,TESTNAME ## LeftConstant,OP)         \
  BINARY_OPRC_TEST(FIXTURENAME,TESTNAME ## RightConstant,OP)

#define RELOP_OP2_TEST(FIXTURENAME,TESTNAME,OP)                     \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {           \
    auto a_dtay = this->a_dtay_;                        \
    auto b_dtay = this->b_dtay_;                        \
    bool r1 = a_dtay OP b_dtay;                         \
    bool r2 = a_dtay.coeff(0) OP b_dtay.coeff(0);       \
    ASSERT_TRUE(r1 == r2);                              \
  }

#define RELOP_OPLC_TEST(FIXTURENAME,TESTNAME,OP)                    \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {           \
    auto b_dtay = this->b_dtay_;                        \
    double val = this->urand.number();                  \
    bool r1 = val OP b_dtay;                            \
    bool r2 = val OP b_dtay.coeff(0);                   \
    ASSERT_TRUE(r1 == r2);                              \
  }

#define RELOP_OPRC_TEST(FIXTURENAME,TESTNAME,OP)                    \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {           \
    auto a_dtay = this->b_dtay_;                        \
    double val = this->urand.number();                  \
    bool r1 = a_dtay OP val;                            \
    bool r2 = a_dtay.coeff(0) OP val;                   \
    ASSERT_TRUE(r1 == r2);                              \
  }

#define RELOP_OP_TEST(FIXTURENAME,TESTNAME,OP)                      \
  RELOP_OP2_TEST(FIXTURENAME,TESTNAME,OP)                           \
  RELOP_OPLC_TEST(FIXTURENAME,TESTNAME ## LeftConstant,OP)          \
  RELOP_OPRC_TEST(FIXTURENAME,TESTNAME ## RightConstant,OP)

#define BINARY_FUNC2_TEST(FIXTURENAME,TESTNAME,FUNC)    \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {   \
    auto a_dtay = this->a_dtay_;            \
    auto b_dtay = this->b_dtay_;            \
    auto c_dtay = this->c_dtay_;            \
    auto d = this->d_;                      \
    auto X = this->X_;                      \
    auto Y = this->Y_;                      \
    c_dtay = FUNC (a_dtay, b_dtay);         \
    trace_on(0);                            \
    adouble aa, ab, ac;                     \
    aa <<= X[0][0];                         \
    ab <<= X[1][0];                         \
    ac = FUNC (aa, ab);                     \
    ac >>= Y[0][0];                         \
    trace_off();                            \
    forward(0,1,2,d,0,X,Y);                 \
    COMPARE_POLYS(c_dtay,Y[0]);             \
  }

#define BINARY_FUNCRC_TEST(FIXTURENAME,TESTNAME,FUNC)   \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {   \
    auto a_dtay = this->a_dtay_;            \
    auto c_dtay = this->c_dtay_;            \
    auto d = this->d_;                      \
    auto X = this->X_;                      \
    auto Y = this->Y_;                      \
    double val = this->urand.number();      \
    c_dtay = FUNC (a_dtay, val);            \
    trace_on(0);                            \
    adouble aa, ac;                         \
    aa <<= X[0][0];                         \
    ac = FUNC (aa, val);                    \
    ac >>= Y[0][0];                         \
    trace_off();                            \
    forward(0,1,1,d,0,X,Y);                 \
    COMPARE_POLYS(c_dtay,Y[0]);             \
  }

#define BINARY_FUNCLC_TEST(FIXTURENAME,TESTNAME,FUNC)   \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {   \
    auto a_dtay = this->a_dtay_;            \
    auto c_dtay = this->c_dtay_;            \
    auto d = this->d_;                      \
    auto X = this->X_;                      \
    auto Y = this->Y_;                      \
    double val = this->urand.number();      \
    c_dtay = FUNC (val, a_dtay);            \
    trace_on(0);                            \
    adouble aa, ac;                         \
    aa <<= X[0][0];                         \
    ac = FUNC (val, aa);                    \
    ac >>= Y[0][0];                         \
    trace_off();                            \
    forward(0,1,1,d,0,X,Y);                 \
    COMPARE_POLYS(c_dtay,Y[0]);             \
  }

#define BINARY_FUNC_TEST(FIXTURENAME,TESTNAME,FUNC)                         \
  BINARY_FUNC2_TEST(FIXTURENAME,TESTNAME,FUNC)                              \
  BINARY_FUNCLC_TEST(FIXTURENAME,TESTNAME ## LeftConstant,FUNC)             \
  BINARY_FUNCRC_TEST(FIXTURENAME,TESTNAME ## RightConstant,FUNC)

#define UNARY_OP_TEST(FIXTURENAME,TESTNAME,OP)                  \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {       \
    auto a_dtay = this->a_dtay_;                    \
    auto c_dtay = this->c_dtay_;                    \
    auto d = this->d_;                              \
    auto X = this->X_;                              \
    auto Y = this->Y_;                              \
    c_dtay = OP a_dtay;                             \
    trace_on(0);                                    \
    adouble aa, ac;                                 \
    aa <<= X[0][0];                                 \
    ac = OP aa;                                     \
    ac >>= Y[0][0];                                 \
    trace_off();                                    \
    forward(0,1,1,d,0,X,Y);                         \
    COMPARE_POLYS(c_dtay,Y[0]);                     \
  }

#define UNARY_FUNC_TEST(FIXTURENAME,TESTNAME,FUNC)              \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {       \
    auto a_dtay = this->a_dtay_;                    \
    auto c_dtay = this->c_dtay_;                    \
    auto d = this->d_;                              \
    auto X = this->X_;                              \
    auto Y = this->Y_;                              \
    c_dtay = FUNC (a_dtay);                         \
    trace_on(0);                                    \
    adouble aa, ac;                                 \
    aa <<= X[0][0];                                 \
    ac = FUNC (aa);                                 \
    ac >>= Y[0][0];                                 \
    trace_off();                                    \
    forward(0,1,1,d,0,X,Y);                         \
    COMPARE_POLYS(c_dtay,Y[0]);                     \
  }

#define UNARY_ASSIGNOP2_TEST(FIXTURENAME,TESTNAME,OP)           \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {       \
    auto a_dtay = this->a_dtay_;                    \
    auto b_dtay = this->b_dtay_;                    \
    auto c_dtay = this->c_dtay_;                    \
    auto d = this->d_;                              \
    auto X = this->X_;                              \
    auto Y = this->Y_;                              \
    c_dtay = a_dtay;                                \
    c_dtay OP b_dtay;                               \
    trace_on(0);                                    \
    adouble aa, ab, ac;                             \
    aa <<= X[0][0];                                 \
    ab <<= X[1][0];                                 \
    ac = aa;                                        \
    ac OP ab;                                       \
    ac >>= Y[0][0];                                 \
    trace_off();                                    \
    forward(0,1,2,d,0,X,Y);                         \
    COMPARE_POLYS(c_dtay,Y[0]);                     \
  }

#define UNARY_ASSIGNOPRC_TEST(FIXTURENAME,TESTNAME,OP)          \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {       \
    auto a_dtay = this->a_dtay_;                    \
    auto c_dtay = this->c_dtay_;                    \
    auto d = this->d_;                              \
    auto X = this->X_;                              \
    auto Y = this->Y_;                              \
    double val = this->urand.number();              \
    c_dtay = a_dtay;                                \
    c_dtay OP val;                                  \
    trace_on(0);                                    \
    adouble aa, ac;                                 \
    aa <<= X[0][0];                                 \
    ac = aa;                                        \
    ac OP val;                                      \
    ac >>= Y[0][0];                                 \
    trace_off();                                    \
    forward(0,1,1,d,0,X,Y);                         \
    COMPARE_POLYS(c_dtay,Y[0]);                     \
  }

#define UNARY_ASSIGNOPLC_TEST(FIXTURENAME,TESTNAME,OP)          \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                         \
    auto a_dtay = this->a_dtay_;                    \
    auto c_dtay = this->c_dtay_;                    \
    auto d = this->d_;                              \
    auto X = this->X_;                              \
    auto Y = this->Y_;                              \
    double val = this->urand.number();              \
    c_dtay = val;                                   \
    c_dtay OP a_dtay;                               \
    trace_on(0);                                    \
    adouble aa, ac;                                 \
    aa <<= X[0][0];                                 \
    ac = val;                                       \
    ac OP aa;                                       \
    ac >>= Y[0][0];                                 \
    trace_off();                                    \
    forward(0,1,1,d,0,X,Y);                         \
    COMPARE_POLYS(c_dtay,Y[0]);                     \
  }

#define UNARY_ASSIGNOP_TEST(FIXTURENAME,TESTNAME,OP)                \
  UNARY_ASSIGNOP2_TEST(FIXTURENAME,TESTNAME,OP)                     \
  UNARY_ASSIGNOPLC_TEST(FIXTURENAME,TESTNAME ## LeftConstant,OP)    \
  UNARY_ASSIGNOPRC_TEST(FIXTURENAME,TESTNAME ## RightConstant,OP)

BINARY_OP_TEST(TaylorOpsUnitTest, testAddition, +)
BINARY_OP_TEST(TaylorOpsUnitTest, testSubtraction, -)
BINARY_OP_TEST(TaylorOpsUnitTest, testMultiplication, *)
BINARY_OP_TEST(TaylorOpsUnitTest, testDivision, /)

RELOP_OP_TEST(TaylorOpsUnitTest, testEquals, ==)
RELOP_OP_TEST(TaylorOpsUnitTest, testNotEquals, !=)
RELOP_OP_TEST(TaylorOpsUnitTest, testLessThanOrEquals, <=)
RELOP_OP_TEST(TaylorOpsUnitTest, testGreaterThanOrEquals, >=)
RELOP_OP_TEST(TaylorOpsUnitTest, testLessThan, <)
RELOP_OP_TEST(TaylorOpsUnitTest, testGreaterThan, >)

BINARY_FUNC_TEST(TaylorOpsUnitTest, testPow, pow)

UNARY_OP_TEST(TaylorOpsUnitTest, testUnaryPlus, +)
UNARY_OP_TEST(TaylorOpsUnitTest, testUnaryMinus, -)

UNARY_FUNC_TEST(TaylorOpsUnitTest, testExp, exp)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testLog, log)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testLog10, log10)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testSqrt, sqrt)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testCos, cos)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testSin, sin)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testTan, tan)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testACos, acos)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testASin, asin)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testATan, atan)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testCosh, cosh)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testSinh, sinh)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testTanh, tanh)
UNARY_FUNC_TEST(TaylorOpsUnitTest, testFAbs, fabs)

UNARY_ASSIGNOP_TEST(TaylorOpsUnitTest, testPlusEquals, +=)
UNARY_ASSIGNOP_TEST(TaylorOpsUnitTest, testMinusEquals, -=)
UNARY_ASSIGNOP_TEST(TaylorOpsUnitTest, testTimesEquals, *=)
UNARY_ASSIGNOP_TEST(TaylorOpsUnitTest, testDivideEquals, /=)

TYPED_TEST_P(TaylorOpsUnitTest, testComposite1) {
  auto a_dtay = this->a_dtay_;
  auto b_dtay = this->b_dtay_;
  auto c_dtay = this->c_dtay_;
  auto d = this->d_;
  auto X = this->X_;
  auto Y = this->Y_;
  c_dtay = this->composite1(a_dtay, b_dtay);
  trace_on(0);
  adouble aa, ab, ac;
  aa <<= X[0][0];
  ab <<= X[1][0];
  ac = this->composite1(aa,ab);
  ac >>= Y[0][0];
  trace_off();
  forward(0,1,2,d,0,X,Y);
  COMPARE_POLYS(c_dtay,Y[0]);
}

TYPED_TEST_P(TaylorOpsUnitTest, testDiff1) {
  typedef decltype(this->a_dtay_) TaylorType;
  auto a_dtay = this->a_dtay_;
  auto d = this->d_;
  TaylorType a_diff1 = diff(a_dtay);
  TaylorType a_diff2(d-1, 0.0);
  for (int i=1; i<=d; ++i)
    a_diff2.fastAccessCoeff(i-1) = a_dtay.fastAccessCoeff(i)*i;
  COMPARE_TAYS(a_diff1, a_diff2);
}

TYPED_TEST_P(TaylorOpsUnitTest, testDiff3) {
  typedef decltype(this->a_dtay_) TaylorType;
  auto a_dtay = this->a_dtay_;
  TaylorType a_diff1 = diff(a_dtay, 3);
  TaylorType a_diff2 = diff( diff( diff(a_dtay) ) );
  COMPARE_TAYS(a_diff1, a_diff2);
}

REGISTER_TYPED_TEST_SUITE_P(
  TaylorOpsUnitTest,
  testAddition,
  testAdditionLeftConstant,
  testAdditionRightConstant,
  testSubtraction,
  testSubtractionLeftConstant,
  testSubtractionRightConstant,
  testMultiplication,
  testMultiplicationLeftConstant,
  testMultiplicationRightConstant,
  testDivision,
  testDivisionLeftConstant,
  testDivisionRightConstant,
  testEquals,
  testEqualsLeftConstant,
  testEqualsRightConstant,
  testNotEquals,
  testNotEqualsLeftConstant,
  testNotEqualsRightConstant,
  testLessThanOrEquals,
  testLessThanOrEqualsLeftConstant,
  testLessThanOrEqualsRightConstant,
  testGreaterThanOrEquals,
  testGreaterThanOrEqualsLeftConstant,
  testGreaterThanOrEqualsRightConstant,
  testLessThan,
  testLessThanLeftConstant,
  testLessThanRightConstant,
  testGreaterThan,
  testGreaterThanLeftConstant,
  testGreaterThanRightConstant,
  testUnaryPlus,
  testUnaryMinus,
  testExp,
  testLog,
  testLog10,
  testSqrt,
  testCos,
  testSin,
  testTan,
  testACos,
  testASin,
  testATan,
  testCosh,
  testSinh,
  testTanh,
  testFAbs,
  testPlusEquals,
  testPlusEqualsLeftConstant,
  testPlusEqualsRightConstant,
  testMinusEquals,
  testMinusEqualsLeftConstant,
  testMinusEqualsRightConstant,
  testTimesEquals,
  testTimesEqualsLeftConstant,
  testTimesEqualsRightConstant,
  testDivideEquals,
  testDivideEqualsLeftConstant,
  testDivideEqualsRightConstant,
  testPow,
  testPowLeftConstant,
  testPowRightConstant,
  testComposite1,
  testDiff1,
  testDiff3);

BINARY_FUNC_TEST(TaylorMaxMinUnitTest, testMax, max)
BINARY_FUNC_TEST(TaylorMaxMinUnitTest, testMin, min)

REGISTER_TYPED_TEST_SUITE_P(
  TaylorMaxMinUnitTest,
  testMax,
  testMaxLeftConstant,
  testMaxRightConstant,
  testMin,
  testMinLeftConstant,
  testMinRightConstant);

#endif // TAYLORUNITTESTS_HPP
