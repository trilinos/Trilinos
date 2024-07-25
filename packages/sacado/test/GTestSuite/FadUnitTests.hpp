// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef FADUNITTESTS_HPP
#define FADUNITTESTS_HPP

// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Random.hpp"

// Fad includes
#include "Fad/fad.h"

// gtest includes
#include <gtest/gtest.h>

#include "GTestUtils.hpp"

// A class for testing each Fad operation
template <typename FadType>
class FadOpsUnitTest : public ::testing::Test {
protected:
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  // DFad variables
  FadType a_dfad_, b_dfad_, c_dfad_;

  // Fad variables
  FAD::Fad<ScalarType> a_fad_, b_fad_, c_fad_;

  // Random number generator
  Sacado::Random<ScalarType> urand;

  // Number of derivative components
  int n;

  // Tolerances to which fad objects should be the same
  ScalarType tol_a, tol_r;

  FadOpsUnitTest() : urand(), n(5), tol_a(1.0e-15), tol_r(1.0e-12) {}

  void SetUp() override {
    ScalarType val;

    val = urand.number();
    a_dfad_ = FadType(n,val);
    a_fad_ = FAD::Fad<ScalarType>(n,val);

    val = urand.number();
    b_dfad_ = FadType(n,val);
    b_fad_ = FAD::Fad<ScalarType>(n,val);

    for (int i=0; i<n; i++) {
      val = urand.number();
      a_dfad_.fastAccessDx(i) = val;
      a_fad_.fastAccessDx(i) = val;

      val = urand.number();
      b_dfad_.fastAccessDx(i) = val;
      b_fad_.fastAccessDx(i) = val;
    }

    val = 0.0;
    c_dfad_ = FadType(n, val);
    c_fad_ = FAD::Fad<ScalarType>(n,val);
  }

  void TearDown() override {}

  template <typename ScalarT>
  ScalarT composite1(const ScalarT& a, const ScalarT& b) {
    ScalarT t1 = 3. * a + sin(b) / log(fabs(a - b * 7.));
    ScalarT t2 = 1.0e3;
    ScalarT t3 = 5.7e4;
    ScalarT t4 = 3.2e5;
    t1 *= cos(a + exp(t1)) / 6. - tan(t1*sqrt(abs(a * log10(abs(b)))));
    t1 -= acos((6.+asin(pow(fabs(a),b)/t2))/t3) * asin(pow(fabs(b),2.)*1.0/t4) * atan((b*pow(2.,log(abs(a))))/(t3*t4));
    t1 /= cosh(b - 0.7) + 7.*sinh(t1 + 0.8)*tanh(9./a) - 9.;
    t1 += pow(abs(a*4.),b-8.)/cos(a*b*a);

    return t1;
  }

}; // class FadOpsUnitTest

TYPED_TEST_SUITE_P(FadOpsUnitTest);

#define BINARY_OP_TEST(FIXTURENAME,TESTNAME,OP)        \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                \
    auto a_dfad = this->a_dfad_;                       \
    auto b_dfad = this->b_dfad_;                       \
    auto c_dfad = this->c_dfad_;                       \
    auto a_fad = this->a_fad_;                         \
    auto b_fad = this->b_fad_;                         \
    auto c_fad = this->c_fad_;                         \
    c_dfad = a_dfad OP b_dfad;                         \
    c_fad = a_fad OP b_fad;                            \
    COMPARE_FADS(c_dfad, c_fad);                       \
                                                       \
    double val = this->urand.number();                 \
    c_dfad = a_dfad OP val;                            \
    c_fad = a_fad OP val;                              \
    COMPARE_FADS(c_dfad, c_fad);                       \
                                                       \
    c_dfad = val OP b_dfad;                            \
    c_fad = val OP b_fad;                              \
    COMPARE_FADS(c_dfad, c_fad);                       \
  }

#define RELOP_TEST(FIXTURENAME,TESTNAME,OP)     \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {         \
    auto a_dfad = this->a_dfad_;                \
    auto b_dfad = this->b_dfad_;                \
    auto a_fad = this->a_fad_;                  \
    auto b_fad = this->b_fad_;                  \
    bool r1 = a_dfad OP b_dfad;                 \
    bool r2 = a_fad OP b_fad;                   \
    ASSERT_TRUE(r1 == r2);                      \
                                                \
    double val = this->urand.number();          \
    r1 = a_dfad OP val;                         \
    r2 = a_fad OP val;                          \
    ASSERT_TRUE(r1 == r2);                      \
                                                \
    r1 = val OP b_dfad;                         \
    r2 = val OP b_fad;                          \
    ASSERT_TRUE(r1 == r2);                      \
  }

#define BINARY_FUNC_TEST(FIXTURENAME,TESTNAME,FUNC)    \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                \
    auto a_dfad = this->a_dfad_;                       \
    auto b_dfad = this->b_dfad_;                       \
    auto c_dfad = this->c_dfad_;                       \
    auto a_fad = this->a_fad_;                         \
    auto b_fad = this->b_fad_;                         \
    auto c_fad = this->c_fad_;                         \
    c_dfad = FUNC (a_dfad,b_dfad);                     \
    c_fad = FUNC (a_fad,b_fad);                        \
    COMPARE_FADS(c_dfad, c_fad);                       \
                                                       \
    double val = this->urand.number();                 \
    c_dfad = FUNC (a_dfad,val);                        \
    c_fad = FUNC (a_fad,val);                          \
    COMPARE_FADS(c_dfad, c_fad);                       \
                                                       \
    c_dfad = FUNC (val,b_dfad);                        \
    c_fad = FUNC (val,b_fad);                          \
    COMPARE_FADS(c_dfad, c_fad);                       \
  }

#define UNARY_OP_TEST(FIXTURENAME,TESTNAME,OP)           \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                  \
    auto a_dfad = this->a_dfad_;                         \
    auto c_dfad = this->c_dfad_;                         \
    auto a_fad = this->a_fad_;                           \
    auto c_fad = this->c_fad_;                           \
    c_dfad = OP a_dfad;                                  \
    c_fad = OP a_fad;                                    \
    COMPARE_FADS(c_dfad, c_fad);                         \
  }

#define UNARY_FUNC_TEST(FIXTURENAME,TESTNAME,FUNC)       \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                  \
    auto a_dfad = this->a_dfad_;                         \
    auto c_dfad = this->c_dfad_;                         \
    auto a_fad = this->a_fad_;                           \
    auto c_fad = this->c_fad_;                           \
    c_dfad = FUNC (a_dfad);                              \
    c_fad = FUNC (a_fad);                                \
    COMPARE_FADS(c_dfad, c_fad);                         \
  }

#define UNARY_ASSIGNOP_TEST(FIXTURENAME,TESTNAME,OP)    \
  TYPED_TEST_P(FIXTURENAME, TESTNAME){                  \
    auto a_dfad = this->a_dfad_;                        \
    auto b_dfad = this->b_dfad_;                        \
    auto c_dfad = this->c_dfad_;                        \
    auto a_fad = this->a_fad_;                          \
    auto b_fad = this->b_fad_;                          \
    auto c_fad = this->c_fad_;                          \
    c_dfad = b_dfad;                                    \
    c_fad = b_fad;                                      \
    c_dfad OP a_dfad;                                   \
    c_fad OP a_fad;                                     \
    COMPARE_FADS(c_dfad, c_fad);                        \
                                                        \
    double val = this->urand.number();                  \
    c_dfad OP val;                                      \
    c_fad OP val;                                       \
    COMPARE_FADS(c_dfad, c_fad);                        \
  }

BINARY_OP_TEST(FadOpsUnitTest, testAddition, +)
BINARY_OP_TEST(FadOpsUnitTest, testSubtraction, -)
BINARY_OP_TEST(FadOpsUnitTest, testMultiplication, *)
BINARY_OP_TEST(FadOpsUnitTest, testDivision, /)

RELOP_TEST(FadOpsUnitTest, testEquals, ==)
RELOP_TEST(FadOpsUnitTest, testNotEquals, !=)
RELOP_TEST(FadOpsUnitTest, testLessThanOrEquals, <=)
RELOP_TEST(FadOpsUnitTest, testGreaterThanOrEquals, >=)
RELOP_TEST(FadOpsUnitTest, testLessThan, <)
RELOP_TEST(FadOpsUnitTest, testGreaterThan, >)

BINARY_FUNC_TEST(FadOpsUnitTest, testPow, pow)

UNARY_OP_TEST(FadOpsUnitTest, testUnaryPlus, +)
UNARY_OP_TEST(FadOpsUnitTest, testUnaryMinus, -)

UNARY_FUNC_TEST(FadOpsUnitTest, testExp, exp)
UNARY_FUNC_TEST(FadOpsUnitTest, testLog, log)
UNARY_FUNC_TEST(FadOpsUnitTest, testLog10, log10)
UNARY_FUNC_TEST(FadOpsUnitTest, testSqrt, sqrt)
UNARY_FUNC_TEST(FadOpsUnitTest, testCos, cos)
UNARY_FUNC_TEST(FadOpsUnitTest, testSin, sin)
UNARY_FUNC_TEST(FadOpsUnitTest, testTan, tan)
UNARY_FUNC_TEST(FadOpsUnitTest, testACos, acos)
UNARY_FUNC_TEST(FadOpsUnitTest, testASin, asin)
UNARY_FUNC_TEST(FadOpsUnitTest, testATan, atan)
UNARY_FUNC_TEST(FadOpsUnitTest, testCosh, cosh)
UNARY_FUNC_TEST(FadOpsUnitTest, testSinh, sinh)
UNARY_FUNC_TEST(FadOpsUnitTest, testTanh, tanh)
UNARY_FUNC_TEST(FadOpsUnitTest, testAbs, abs)
UNARY_FUNC_TEST(FadOpsUnitTest, testFAbs, fabs)

UNARY_ASSIGNOP_TEST(FadOpsUnitTest, testPlusEquals, +=)
UNARY_ASSIGNOP_TEST(FadOpsUnitTest, testMinusEquals, -=)
UNARY_ASSIGNOP_TEST(FadOpsUnitTest, testTimesEquals, *=)
UNARY_ASSIGNOP_TEST(FadOpsUnitTest, testDivideEquals, /=)

TYPED_TEST_P(FadOpsUnitTest, testMax) {
  typedef decltype(this->a_dfad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  ScalarType val;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto c_dfad = this->c_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_dfad = a_dfad + 1.0;
  c_dfad = max(aa_dfad, a_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad, aa_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad+1.0, a_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad, a_dfad+1.0);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  val = a_dfad.val() + 1;
  c_dfad = max(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);

  val = a_dfad.val() - 1;
  c_dfad = max(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), a_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), a_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() + 1;
  c_dfad = max(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);

  val = b_dfad.val() - 1;
  c_dfad = max(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), b_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), b_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}

TYPED_TEST_P(FadOpsUnitTest, testMin) {
  typedef decltype(this->a_dfad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  ScalarType val;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto c_dfad = this->c_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_dfad = a_dfad - 1.0;
  c_dfad = min(aa_dfad, a_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = min(a_dfad, aa_dfad);
  COMPARE_VALUES(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  val = a_dfad.val() - 1;
  c_dfad = min(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);

  val = a_dfad.val() + 1;
  c_dfad = min(a_dfad, val);
  COMPARE_VALUES(c_dfad.val(), a_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), a_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() - 1;
  c_dfad = min(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), val);
  for (int i=0; i<this->n; i++)
    COMPARE_VALUES(c_dfad.dx(i), 0.0);

  val = b_dfad.val() + 1;
  c_dfad = min(val, b_dfad);
  COMPARE_VALUES(c_dfad.val(), b_dfad.val());
  for (int i=0; i<this->n; i++) {
    COMPARE_VALUES(c_dfad.dx(i), b_dfad.dx(i));
    COMPARE_VALUES(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}

TYPED_TEST_P(FadOpsUnitTest, testComposite1) {
  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto c_dfad = this->c_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  c_dfad = this->composite1(a_dfad, b_dfad);
  c_fad = this->composite1(a_fad, b_fad);
  COMPARE_FADS(c_dfad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest, testPlusLR) {
  typedef decltype(this->a_dfad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadType aa_dfad = a_dfad;
  FAD::Fad<ScalarType> aa_fad = a_fad;
  aa_dfad = 1.0;
  aa_fad = 1.0;
  aa_dfad = aa_dfad + b_dfad;
  aa_fad = aa_fad + b_fad;
  COMPARE_FADS(aa_dfad, aa_fad);
}

TYPED_TEST_P(FadOpsUnitTest, testMinusLR) {
  typedef decltype(this->a_dfad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadType aa_dfad = a_dfad;
  FAD::Fad<ScalarType> aa_fad = a_fad;
  aa_dfad = 1.0;
  aa_fad = 1.0;
  aa_dfad = aa_dfad - b_dfad;
  aa_fad = aa_fad - b_fad;
  COMPARE_FADS(aa_dfad, aa_fad);
}

TYPED_TEST_P(FadOpsUnitTest, testTimesLR) {
  typedef decltype(this->a_dfad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadType aa_dfad = a_dfad;
  FAD::Fad<ScalarType> aa_fad = a_fad;
  aa_dfad = 2.0;
  aa_fad = 2.0;
  aa_dfad = aa_dfad * b_dfad;
  aa_fad = aa_fad * b_fad;
  COMPARE_FADS(aa_dfad, aa_fad);
}

TYPED_TEST_P(FadOpsUnitTest, testDivideLR) {
  typedef decltype(this->a_dfad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadType aa_dfad = a_dfad;
  FAD::Fad<ScalarType> aa_fad = a_fad;
  aa_dfad = 2.0;
  aa_fad = 2.0;
  aa_dfad = aa_dfad / b_dfad;
  aa_fad = aa_fad / b_fad;
  COMPARE_FADS(aa_dfad, aa_fad);
}

  // Check various corner cases for pow()
TYPED_TEST_P(FadOpsUnitTest, testPowConstB) {
  typedef decltype(this->a_dfad_) FadType;

  FadType a, b, c, cc;

  // Constant b
  a = FadType(this->n,1.2345);
  for (int i=0; i<this->n; ++i)
    a.fastAccessDx(i) = this->urand.number();
  b = 3.456;
  c = pow(a, b);
  cc = FadType(this->n, pow(a.val(),b.val()));
  for (int i=0; i<this->n; ++i)
    cc.fastAccessDx(i) = b.val()*pow(a.val(),b.val()-1)*a.dx(i);
  COMPARE_FADS(c, cc);

  // Constant scalar b
  c = pow(a, b.val());
  COMPARE_FADS(c, cc);

  // Constant b == 0
  b = 0.0;
  c = pow(a, b);
  cc.val() = 1.0;
  for (int i=0; i<this->n; ++i)
    cc.fastAccessDx(i) = 0.0;
  COMPARE_FADS(c, cc);

  // Constant scalar b == 0
  c = pow(a, b.val());
  COMPARE_FADS(c, cc);

  // a == 0 and constant b
  a.val() = 0.0;
  b = 3.456;
  c = pow(a, b);
  cc.val() = 0.0;
  for (int i=0; i<this->n; ++i)
    cc.fastAccessDx(i) = 0.0;
  COMPARE_FADS(c, cc);

  // a == 0 and constant scalar b
  c = pow(a, b.val());
  COMPARE_FADS(c, cc);

  // a == 0 and b == 0
  b = 0.0;
  c = pow(a, b);
  cc.val() = 1.0;
  for (int i=0; i<this->n; ++i)
    cc.fastAccessDx(i) = 0.0;
  COMPARE_FADS(c, cc);

  // a == 0 and scalar b == 0
  c = pow(a, b.val());
  COMPARE_FADS(c, cc);

  // a == 0 and b == 1
  b = 1.0;
  c = pow(a, b);
  cc = a;
  if (!Sacado::IsStaticallySized<FadType>::value) {
    COMPARE_FADS(c, cc);
  }
  c = pow(a, b.val());
  COMPARE_FADS(c, cc);

  // a == 0 and b == 2
  b = 2.0;
  c = pow(a, b);
  cc = a*a;
  if (!Sacado::IsStaticallySized<FadType>::value) {
    COMPARE_FADS(c, cc);
  }
  c = pow(a, b.val());
  COMPARE_FADS(c, cc);
}

REGISTER_TYPED_TEST_SUITE_P(
  FadOpsUnitTest,
  testAddition,
  testSubtraction,
  testMultiplication,
  testDivision,
  testEquals,
  testNotEquals,
  testLessThanOrEquals,
  testGreaterThanOrEquals,
  testLessThan,
  testGreaterThan,
  testPow,
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
  testAbs,
  testFAbs,
  testPlusEquals,
  testMinusEquals,
  testTimesEquals,
  testDivideEquals,
  testMax,
  testMin,
  testComposite1,
  testPlusLR,
  testMinusLR,
  testTimesLR,
  testDivideLR,
  testPowConstB);

#endif // FADUNITTESTS_HPP
