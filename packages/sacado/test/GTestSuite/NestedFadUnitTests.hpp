// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef NESTED_FADUNITTESTS_HPP
#define NESTED_FADUNITTESTS_HPP

// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Random.hpp"

// Fad includes
#include "Fad/fad.h"

// gtest includes
#include <gtest/gtest.h>

#include "GTestUtils.hpp"

// A class for testing each Fad operation
template <class FadFadType>
class FadFadOpsUnitTest : public ::testing::Test {
protected:
  typedef typename Sacado::ScalarType<FadFadType>::type ScalarType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;

  // DFad variables
  FadFadType a_dfad_, b_dfad_, c_dfad_;

  // Fad variables
  FAD::Fad<FadType> a_fad_, b_fad_, c_fad_;

  // Random number generator
  Sacado::Random<ScalarType> urand;

  // Number of derivative components
  int n1, n2;

  // Tolerances to which fad objects should be the same
  ScalarType tol_a, tol_r;

  // Set the random number generator with a specific seed to prevent NaN's
  // in the second derivatives, likely due to overflow
  FadFadOpsUnitTest() :
    urand(0.0, 1.0, 123456), n1(5), n2(3), tol_a(1.0e-15), tol_r(1.0e-14) {}

  void SetUp() override {
    ScalarType val;

    val = urand.number();
    a_dfad_ = FadFadType(n1,FadType(n2,val));
    a_fad_ = FAD::Fad<FadType>(n1,FadType(n2,val));

    val = urand.number();
    b_dfad_ = FadFadType(n1,FadType(n2,val));
    b_fad_ = FAD::Fad<FadType>(n1,FadType(n2,val));

    for (int j=0; j<n2; j++) {
      ScalarType val2;
      val2 = urand.number();
      a_dfad_.val().fastAccessDx(j) = val2;
      a_fad_.val().fastAccessDx(j) = val2;

      val2 = urand.number();
      b_dfad_.val().fastAccessDx(j) = val2;
      b_fad_.val().fastAccessDx(j) = val2;
    }

    for (int i=0; i<n1; i++) {
      val = urand.number();
      a_dfad_.fastAccessDx(i) = FadType(n2,val);
      a_fad_.fastAccessDx(i) = FadType(n2,val);

      val = urand.number();
      b_dfad_.fastAccessDx(i) = FadType(n2,val);
      b_fad_.fastAccessDx(i) = FadType(n2,val);

      for (int j=0; j<n2; j++) {
        ScalarType val2;
        val2 = urand.number();
        a_dfad_.fastAccessDx(i).fastAccessDx(j) = val2;
        a_fad_.fastAccessDx(i).fastAccessDx(j) = val2;

        val2 = urand.number();
        b_dfad_.fastAccessDx(i).fastAccessDx(j) = val2;
        b_fad_.fastAccessDx(i).fastAccessDx(j) = val2;
      }
    }

    val = 0.0;
    c_dfad_ = FadFadType(n1, FadType(n2,val));
    c_fad_ = FAD::Fad<FadType>(n1,FadType(n2,val));
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

  template <typename ScalarT>
  ScalarT composite1_fad(const ScalarT& a, const ScalarT& b) {
    ScalarT t1 = FadType(3.) * a + sin(b) / log(fabs(a - b * FadType(7.)));
    ScalarT t2 = FadType(1.0e3);
    ScalarT t3 = FadType(7e4);
    ScalarT t4 = FadType(3.2e5);
    t1 *= cos(a + exp(t1)) / FadType(6.) - tan(t1*sqrt(abs(a * log10(abs(b)))));
    t1 -= acos((FadType(6.)+asin(pow(fabs(a),b)/t2))/t3) * asin(pow(fabs(b),FadType(2.))*FadType(1.0)/t4) * atan((b*pow(FadType(2.),log(abs(a))))/(t3*t4));
    t1 /= cosh(b - FadType(0.7)) + FadType(7.)*sinh(t1 + FadType(0.8))*tanh(FadType(9.)/a) - FadType(9.);
    t1 += pow(abs(a*FadType(4.)),b-FadType(8.))/cos(a*b*a);

    return t1;
  }

}; // class FadFadOpsUnitTest

TYPED_TEST_SUITE_P(FadFadOpsUnitTest);

#define BINARY_OP_TEST(FIXTURENAME,TESTNAME,OP)        \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                \
    typedef decltype(this->a_dfad_) FadFadType;                   \
    typedef typename Sacado::ValueType<FadFadType>::type FadType; \
    auto a_dfad = this->a_dfad_;                       \
    auto b_dfad = this->b_dfad_;                       \
    auto c_dfad = this->c_dfad_;                       \
    auto a_fad = this->a_fad_;                         \
    auto b_fad = this->b_fad_;                         \
    auto c_fad = this->c_fad_;                         \
    c_dfad = a_dfad OP b_dfad;                         \
    c_fad = a_fad OP b_fad;                            \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                \
                                                       \
    double val = this->urand.number();                 \
    c_dfad = a_dfad OP val;                            \
    c_fad = a_fad OP FadType(val);                     \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                \
                                                       \
    c_dfad = val OP b_dfad;                            \
    c_fad = FadType(val) OP b_fad;                     \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                \
  }

#define RELOP_TEST(FIXTURENAME,TESTNAME,OP)     \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {         \
    typedef decltype(this->a_dfad_) FadFadType;                   \
    typedef typename Sacado::ValueType<FadFadType>::type FadType; \
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
    r2 = a_fad OP FadType(val);                 \
    ASSERT_TRUE(r1 == r2);                      \
                                                \
    r1 = val OP b_dfad;                         \
    r2 = FadType(val) OP b_fad;                 \
    ASSERT_TRUE(r1 == r2);                      \
  }

#define BINARY_FUNC_TEST(FIXTURENAME,TESTNAME,FUNC)    \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                \
    typedef decltype(this->a_dfad_) FadFadType;                   \
    typedef typename Sacado::ValueType<FadFadType>::type FadType; \
    auto a_dfad = this->a_dfad_;                       \
    auto b_dfad = this->b_dfad_;                       \
    auto c_dfad = this->c_dfad_;                       \
    auto a_fad = this->a_fad_;                         \
    auto b_fad = this->b_fad_;                         \
    auto c_fad = this->c_fad_;                         \
    c_dfad = FUNC (a_dfad,b_dfad);                     \
    c_fad = FUNC (a_fad,b_fad);                        \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                \
                                                       \
    double val = this->urand.number();                 \
    c_dfad = FUNC (a_dfad,val);                        \
    c_fad = FUNC (a_fad,FadType(val));                 \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                \
                                                       \
    c_dfad = FUNC (val,b_dfad);                        \
    c_fad = FUNC (FadType(val),b_fad);                 \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                \
  }

#define UNARY_OP_TEST(FIXTURENAME,TESTNAME,OP)           \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                  \
    auto a_dfad = this->a_dfad_;                         \
    auto c_dfad = this->c_dfad_;                         \
    auto a_fad = this->a_fad_;                           \
    auto c_fad = this->c_fad_;                           \
    c_dfad = OP a_dfad;                                  \
    c_fad = OP a_fad;                                    \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                  \
  }

#define UNARY_FUNC_TEST(FIXTURENAME,TESTNAME,FUNC)       \
  TYPED_TEST_P(FIXTURENAME, TESTNAME) {                  \
    auto a_dfad = this->a_dfad_;                         \
    auto c_dfad = this->c_dfad_;                         \
    auto a_fad = this->a_fad_;                           \
    auto c_fad = this->c_fad_;                           \
    c_dfad = FUNC (a_dfad);                              \
    c_fad = FUNC (a_fad);                                \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                  \
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
    COMPARE_NESTED_FADS(c_dfad, c_fad);                 \
                                                        \
    double val = this->urand.number();                  \
    c_dfad OP val;                                      \
    c_fad OP val;                                       \
    COMPARE_NESTED_FADS(c_dfad, c_fad);                 \
  }

BINARY_OP_TEST(FadFadOpsUnitTest, testAddition, +)
BINARY_OP_TEST(FadFadOpsUnitTest, testSubtraction, -)
BINARY_OP_TEST(FadFadOpsUnitTest, testMultiplication, *)
BINARY_OP_TEST(FadFadOpsUnitTest, testDivision, /)

RELOP_TEST(FadFadOpsUnitTest, testEquals, ==)
RELOP_TEST(FadFadOpsUnitTest, testNotEquals, !=)
RELOP_TEST(FadFadOpsUnitTest, testLessThanOrEquals, <=)
RELOP_TEST(FadFadOpsUnitTest, testGreaterThanOrEquals, >=)
RELOP_TEST(FadFadOpsUnitTest, testLessThan, <)
RELOP_TEST(FadFadOpsUnitTest, testGreaterThan, >)

BINARY_FUNC_TEST(FadFadOpsUnitTest, testPow, pow)

UNARY_OP_TEST(FadFadOpsUnitTest, testUnaryPlus, +)
UNARY_OP_TEST(FadFadOpsUnitTest, testUnaryMinus, -)

UNARY_FUNC_TEST(FadFadOpsUnitTest, testExp, exp)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testLog, log)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testLog10, log10)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testSqrt, sqrt)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testCos, cos)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testSin, sin)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testTan, tan)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testACos, acos)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testASin, asin)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testATan, atan)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testCosh, cosh)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testSinh, sinh)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testTanh, tanh)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testAbs, abs)
UNARY_FUNC_TEST(FadFadOpsUnitTest, testFAbs, fabs)

UNARY_ASSIGNOP_TEST(FadFadOpsUnitTest, testPlusEquals, +=)
UNARY_ASSIGNOP_TEST(FadFadOpsUnitTest, testMinusEquals, -=)
UNARY_ASSIGNOP_TEST(FadFadOpsUnitTest, testTimesEquals, *=)
UNARY_ASSIGNOP_TEST(FadFadOpsUnitTest, testDivideEquals, /=)

TYPED_TEST_P(FadFadOpsUnitTest, testMax) {
  typedef decltype(this->a_dfad_) FadFadType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;

  FadType val;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto c_dfad = this->c_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadFadType aa_dfad = a_dfad + 1.0;
  c_dfad = max(aa_dfad, a_dfad);
  COMPARE_FADS(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad, aa_dfad);
  COMPARE_FADS(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad+1.0, a_dfad);
  COMPARE_FADS(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad, a_dfad+1.0);
  COMPARE_FADS(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  val = a_dfad.val() + 1;
  c_dfad = max(a_dfad, val);
  COMPARE_FADS(c_dfad.val(), val);
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), FadType(0.0));
  }

  val = a_dfad.val() - 1;
  c_dfad = max(a_dfad, val);
  COMPARE_FADS(c_dfad.val(), a_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), a_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() + 1;
  c_dfad = max(val, b_dfad);
  COMPARE_FADS(c_dfad.val(), val);
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), FadType(0.0));
  }

  val = b_dfad.val() - 1;
  c_dfad = max(val, b_dfad);
  COMPARE_FADS(c_dfad.val(), b_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), b_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}

TYPED_TEST_P(FadFadOpsUnitTest, testMin) {
  typedef decltype(this->a_dfad_) FadFadType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;

  FadType val;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto c_dfad = this->c_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadFadType aa_dfad = a_dfad - 1.0;
  c_dfad = min(aa_dfad, a_dfad);
  COMPARE_FADS(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = min(a_dfad, aa_dfad);
  COMPARE_FADS(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), aa_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  val = a_dfad.val() - 1;
  c_dfad = min(a_dfad, val);
  COMPARE_FADS(c_dfad.val(), val);
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), FadType(0.0));
  }

  val = a_dfad.val() + 1;
  c_dfad = min(a_dfad, val);
  COMPARE_FADS(c_dfad.val(), a_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), a_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() - 1;
  c_dfad = min(val, b_dfad);
  COMPARE_FADS(c_dfad.val(), val);
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), FadType(0.0));
  }

  val = b_dfad.val() + 1;
  c_dfad = min(val, b_dfad);
  COMPARE_FADS(c_dfad.val(), b_dfad.val());
  for (int i=0; i<this->n1; i++) {
    COMPARE_FADS(c_dfad.dx(i), b_dfad.dx(i));
    COMPARE_FADS(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}

TYPED_TEST_P(FadFadOpsUnitTest, testComposite1) {
  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto c_dfad = this->c_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  c_dfad = this->composite1(a_dfad, b_dfad);
  c_fad = this->composite1_fad(a_fad, b_fad);
  COMPARE_NESTED_FADS(c_dfad, c_fad);
}

TYPED_TEST_P(FadFadOpsUnitTest, testPlusLR) {
  typedef decltype(this->a_dfad_) FadFadType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadFadType aa_dfad = a_dfad;
  FAD::Fad< FadType > aa_fad = a_fad;
  aa_dfad = 1.0;
  aa_fad = FadType(1.0);
  aa_dfad = aa_dfad + b_dfad;
  aa_fad = aa_fad + b_fad;
  COMPARE_NESTED_FADS(aa_dfad, aa_fad);
}

TYPED_TEST_P(FadFadOpsUnitTest, testMinusLR) {
  typedef decltype(this->a_dfad_) FadFadType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadFadType aa_dfad = a_dfad;
  FAD::Fad< FadType > aa_fad = a_fad;
  aa_dfad = 1.0;
  aa_fad = FadType(1.0);
  aa_dfad = aa_dfad - b_dfad;
  aa_fad = aa_fad - b_fad;
  COMPARE_NESTED_FADS(aa_dfad, aa_fad);
}

TYPED_TEST_P(FadFadOpsUnitTest, testTimesLR) {
  typedef decltype(this->a_dfad_) FadFadType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadFadType aa_dfad = a_dfad;
  FAD::Fad< FadType > aa_fad = a_fad;
  aa_dfad = 2.0;
  aa_fad = FadType(2.0);
  aa_dfad = aa_dfad * b_dfad;
  aa_fad = aa_fad * b_fad;
  COMPARE_NESTED_FADS(aa_dfad, aa_fad);
}

TYPED_TEST_P(FadFadOpsUnitTest, testDivideLR) {
  typedef decltype(this->a_dfad_) FadFadType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;

  auto a_dfad = this->a_dfad_;
  auto b_dfad = this->b_dfad_;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  FadFadType aa_dfad = a_dfad;
  FAD::Fad< FadType > aa_fad = a_fad;
  aa_dfad = 2.0;
  aa_fad = FadType(2.0);
  aa_dfad = aa_dfad / b_dfad;
  aa_fad = aa_fad / b_fad;
  COMPARE_NESTED_FADS(aa_dfad, aa_fad);
}

  // Check various corner cases for pow()
TYPED_TEST_P(FadFadOpsUnitTest, testPowConstB) {
  typedef decltype(this->a_dfad_) FadFadType;
  typedef typename Sacado::ValueType<FadFadType>::type FadType;
  typedef typename Sacado::ScalarType<FadFadType>::type ScalarType;

  auto a_dfad = this->a_dfad_;

  FadFadType a, b, c, cc;

  // Constant b
  a = a_dfad;
  b = 3.456;
  c = pow(a, b);
  ScalarType f = pow(a.val().val(), b.val().val());
  ScalarType fp = b.val().val()*pow(a.val().val(),b.val().val()-1);
  ScalarType fpp = b.val().val()*(b.val().val()-1)*pow(a.val().val(),b.val().val()-2);
  cc = FadFadType(this->n1,FadType(this->n2,f));
  for (int i=0; i<this->n2; ++i)
    cc.val().fastAccessDx(i) = fp*a.val().dx(i);
  for (int i=0; i<this->n1; ++i) {
    cc.fastAccessDx(i) = FadType(this->n2,fp*a.dx(i).val());
    for (int j=0; j<this->n2; ++j)
      cc.fastAccessDx(i).fastAccessDx(j) = fpp*a.dx(i).val()*a.val().dx(j) + fp*a.dx(i).dx(j);
  }
  COMPARE_NESTED_FADS(c, cc);

  // Constant scalar b
  c = pow(a, b.val());
  COMPARE_NESTED_FADS(c, cc);
  c = pow(a, b.val().val());
  COMPARE_NESTED_FADS(c, cc);

  // Constant b == 0
  b = 0.0;
  c = pow(a, b);
  cc.val() = FadType(this->n2,1.0);
  for (int i=0; i<this->n1; ++i)
    cc.fastAccessDx(i) = 0.0;
  COMPARE_NESTED_FADS(c, cc);

  // Constant scalar b == 0
  c = pow(a, b.val());
  cc.val() = FadType(this->n2,1.0);
  for (int i=0; i<this->n1; ++i)
    cc.fastAccessDx(i) = 0.0;
  COMPARE_NESTED_FADS(c, cc);
  c = pow(a, b.val().val());
  COMPARE_NESTED_FADS(c, cc);

  // a == 0 and constant b as a Fad
  // This only works for DFad/SLFad, because there is no such thing as a
  // constant SFad.
  if (!Sacado::IsStaticallySized<FadType>::value) {
    a.val() = 0.0;
    b = 3.456;
    c = pow(a, b);
    cc.val() = 0.0;
    for (int i=0; i<this->n1; ++i)
      cc.fastAccessDx(i) = FadType(this->n2,0.0);
    COMPARE_NESTED_FADS(c, cc);
  }

  // a == 0 and constant scalar b
  a.val() = 0.0;
  b = 3.456;
  c = pow(a, b.val());
  cc.val() = 0.0;
  for (int i=0; i<this->n1; ++i)
    cc.fastAccessDx(i) = FadType(this->n2,0.0);
  COMPARE_NESTED_FADS(c, cc);
  c = pow(a, b.val().val());
  COMPARE_NESTED_FADS(c, cc);

  // a == 0 and b == 0
  b = 0.0;
  cc.val() = 1.0;
  for (int i=0; i<this->n1; ++i)
    cc.fastAccessDx(i) = 0.0;
  if (!Sacado::IsStaticallySized<FadType>::value) {
    c = pow(a, b);
    COMPARE_NESTED_FADS(c, cc);
    c = pow(a, b.val());
    COMPARE_NESTED_FADS(c, cc);
  }
  c = pow(a, b.val().val());
  COMPARE_NESTED_FADS(c, cc);

  // a nonzero and b == 2
  a = a_dfad;
  a.val().val() = 0.0;
  b = 2.0;
  c = pow(a, b);
  f = pow(a.val().val(), b.val().val());
  fp = b.val().val()*pow(a.val().val(),b.val().val()-1);
  fpp = b.val().val()*(b.val().val()-1)*pow(a.val().val(),b.val().val()-2);
  cc = FadFadType(this->n1,FadType(this->n2,f));
  for (int i=0; i<this->n2; ++i)
    cc.val().fastAccessDx(i) = fp*a.val().dx(i);
  for (int i=0; i<this->n1; ++i) {
    cc.fastAccessDx(i) = FadType(this->n2,fp*a.dx(i).val());
    for (int j=0; j<this->n2; ++j)
      cc.fastAccessDx(i).fastAccessDx(j) = fpp*a.dx(i).val()*a.val().dx(j) + fp*a.dx(i).dx(j);
  }

  // a.val().val() == 0 with a.val().dx() != 0
  if (!Sacado::IsStaticallySized<FadType>::value) {
    COMPARE_NESTED_FADS(c, cc);
    c = pow(a, b.val());
    COMPARE_NESTED_FADS(c, cc);
  }
  c = pow(a, b.val().val());
  COMPARE_NESTED_FADS(c, cc);

  // a.val().val() == 0 and b == 1
  b = 1.0;
  c = pow(a, b);
  cc = a;
  if (!Sacado::IsStaticallySized<FadType>::value) {
    COMPARE_NESTED_FADS(c, cc);
    c = pow(a, b.val());
    COMPARE_NESTED_FADS(c, cc);
  }
  c = pow(a, b.val().val());
  COMPARE_NESTED_FADS(c, cc);

  // a.val().val() == 0 and b == 0
  b = 0.0;
  c = pow(a, b);
  cc.val() = FadType(this->n2, 1.0);
  for (int i=0; i<this->n1; ++i)
    cc.fastAccessDx(i) = 0.0;
  if (!Sacado::IsStaticallySized<FadType>::value) {
    COMPARE_NESTED_FADS(c, cc);
    c = pow(a, b.val());
    COMPARE_NESTED_FADS(c, cc);
  }
  c = pow(a, b.val().val());
  COMPARE_NESTED_FADS(c, cc);

  // a.val().val() == 0, b == 2, nested expression
  c = pow(2.0*a,2.0);
  cc = 4.0*a*a;
  COMPARE_NESTED_FADS(c, cc);
}

REGISTER_TYPED_TEST_SUITE_P(
  FadFadOpsUnitTest,
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

#endif // NESETD_FADUNITTESTS_HPP
