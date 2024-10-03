// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef FADUNITTESTS2_HPP
#define FADUNITTESTS2_HPP

// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Random.hpp"

// gtest includes
#include <gtest/gtest.h>

#include "GTestUtils.hpp"

// A class for testing each Fad operation
template <typename FadType>
class FadOpsUnitTest2 : public ::testing::Test {
protected:
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;

  // DFad variables
  FadType a_fad_, b_fad_, c_fad_;

  // Random number generator
  Sacado::Random<ScalarType> urand;

  // Number of derivative components
  int n_;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

  FadOpsUnitTest2() : urand(), n_(5), tol_a(1.0e-15), tol_r(1.0e-12) {}

  void SetUp() override {
    ScalarType val;

    val = urand.number();
    a_fad_ = FadType(n_,val);

    val = urand.number();
    b_fad_ = FadType(n_,val);

    for (int i=0; i<n_; i++) {
      val = urand.number();
      a_fad_.fastAccessDx(i) = val;

      val = urand.number();
      b_fad_.fastAccessDx(i) = val;
    }

    val = 0.0;
    c_fad_ = FadType(n_, val);
  }

  void TearDown() override {}

}; // class FadOpsUnitTest2

// A class for testing each real Fad operation
// This class tests additional functions that aren't define for complex
// types
template <typename FadType>
class RealFadOpsUnitTest2 : public FadOpsUnitTest2<FadType> {
protected:

  RealFadOpsUnitTest2() : FadOpsUnitTest2<FadType>()  {}

}; // class RealFadOpsUnitTest2

TYPED_TEST_SUITE_P(FadOpsUnitTest2);
TYPED_TEST_SUITE_P(RealFadOpsUnitTest2);

TYPED_TEST_P(FadOpsUnitTest2, testAddition) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = a_fad + b_fad;
  FadType t1(n, a_fad.val()+b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i) + b_fad.dx(i);
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  c_fad = a_fad + val;
  FadType t2(n, a_fad.val()+val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i);
  COMPARE_FADS(c_fad, t2);

  c_fad = val + b_fad;
  FadType t3(n, val+b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);
}

TYPED_TEST_P(FadOpsUnitTest2, testSubtraction) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = a_fad - b_fad;
  FadType t1(n, a_fad.val()-b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i) - b_fad.dx(i);
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  c_fad = a_fad - val;
  FadType t2(n, a_fad.val()-val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i);
  COMPARE_FADS(c_fad, t2);

  c_fad = val - b_fad;
  FadType t3(n, val-b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = -b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);
}

TYPED_TEST_P(FadOpsUnitTest2, testMultiplication) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = a_fad * b_fad;
  FadType t1(n, a_fad.val()*b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)*b_fad.val() + a_fad.val()*b_fad.dx(i);
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  c_fad = a_fad * val;
  FadType t2(n, a_fad.val()*val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i)*val;
  COMPARE_FADS(c_fad, t2);

  c_fad = val * b_fad;
  FadType t3(n, val*b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = val*b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);
}

TYPED_TEST_P(FadOpsUnitTest2, testDivision) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = a_fad / b_fad;
  FadType t1(n, a_fad.val()/b_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) =
      (a_fad.dx(i)*b_fad.val() - a_fad.val()*b_fad.dx(i)) /
      (b_fad.val()*b_fad.val());
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  c_fad = a_fad / val;
  FadType t2(n, a_fad.val()/val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = a_fad.dx(i)/val;
  COMPARE_FADS(c_fad, t2);

  c_fad = val / b_fad;
  FadType t3(n, val/b_fad.val());
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = -val*b_fad.dx(i)/(b_fad.val()*b_fad.val());
  COMPARE_FADS(c_fad, t3);
}

TYPED_TEST_P(FadOpsUnitTest2, testEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  bool r1 = a_fad == b_fad;
  bool r2 = a_fad.val() == b_fad.val();
  ASSERT_TRUE(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = a_fad == val;
  r2 = a_fad.val() == val;
  ASSERT_TRUE(r1 == r2);

  r1 = val == b_fad;
  r2 = val == b_fad.val();
  ASSERT_TRUE(r1 == r2);
}

TYPED_TEST_P(FadOpsUnitTest2, testNotEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  bool r1 = a_fad != b_fad;
  bool r2 = a_fad.val() != b_fad.val();
  ASSERT_TRUE(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = a_fad != val;
  r2 = a_fad.val() != val;
  ASSERT_TRUE(r1 == r2);

  r1 = val != b_fad;
  r2 = val != b_fad.val();
  ASSERT_TRUE(r1 == r2);
}

TYPED_TEST_P(FadOpsUnitTest2, testUnaryPlus) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = +(a_fad);
  FadType t1(n, a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testUnaryMinus) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = -(a_fad);
  FadType t1(n, -a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = -a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testExp) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::exp(a_fad);
  FadType t1(n, std::exp(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::exp(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testLog) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::log(a_fad);
  FadType t1(n, std::log(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/a_fad.val();
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testLog10) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::log10(a_fad);
  FadType t1(n, std::log10(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/(a_fad.val()*std::log(10));
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testSqrt) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::sqrt(a_fad);
  FadType t1(n, std::sqrt(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/(2.*std::sqrt(a_fad.val()));
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testCos) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::cos(a_fad);
  FadType t1(n, std::cos(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = -std::sin(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testSin) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::sin(a_fad);
  FadType t1(n, std::sin(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::cos(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testTan) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::tan(a_fad);
  FadType t1(n, std::tan(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) =
      a_fad.dx(i)/(std::cos(a_fad.val())*std::cos(a_fad.val()));;
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testCosh) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::cosh(a_fad);
  FadType t1(n, std::cosh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::sinh(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testSinh) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::sinh(a_fad);
  FadType t1(n, std::sinh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = std::cosh(a_fad.val())*a_fad.dx(i);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testTanh) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::tanh(a_fad);
  FadType t1(n, std::tanh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) =
      a_fad.dx(i)/(std::cosh(a_fad.val())*std::cosh(a_fad.val()));
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(FadOpsUnitTest2, testPlusEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  FadType t1(n, c_fad.val()+a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = c_fad.dx(i) + a_fad.dx(i);
  c_fad += a_fad;
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  FadType t2(n, c_fad.val()+val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i);
  c_fad += val;
  COMPARE_FADS(c_fad, t2);
}

TYPED_TEST_P(FadOpsUnitTest2, testMinusEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  FadType t1(n, c_fad.val()-a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = c_fad.dx(i) - a_fad.dx(i);
  c_fad -= a_fad;
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  FadType t2(n, c_fad.val()-val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i);
  c_fad -= val;
  COMPARE_FADS(c_fad, t2);
}

TYPED_TEST_P(FadOpsUnitTest2, testTimesEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  FadType t1(n, c_fad.val()*a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = c_fad.dx(i)*a_fad.val() + a_fad.dx(i)*c_fad.val();
  c_fad *= a_fad;
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  FadType t2(n, c_fad.val()*val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i)*val;
  c_fad *= val;
  COMPARE_FADS(c_fad, t2);
}

TYPED_TEST_P(FadOpsUnitTest2, testDivideEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  FadType t1(n, c_fad.val()/a_fad.val());
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) =
      (a_fad.dx(i)*c_fad.val() - c_fad.dx(i)*a_fad.val()) /
      (a_fad.val()*a_fad.val());
  c_fad /= a_fad;
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  FadType t2(n, c_fad.val()/val);
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = c_fad.dx(i)/val;
  c_fad /= val;
  COMPARE_FADS(c_fad, t2);
}

TYPED_TEST_P(FadOpsUnitTest2, testPow) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::pow(a_fad, b_fad);
  FadType t1(n, std::pow(a_fad.val(),b_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) =
      std::pow(a_fad.val(),b_fad.val())*(b_fad.val()*a_fad.dx(i)/a_fad.val() +
                                         std::log(a_fad.val())*b_fad.dx(i));
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  c_fad = std::pow(a_fad, val);
  FadType t2(n, std::pow(a_fad.val(), val));
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) =
      std::pow(a_fad.val(), val)*(val*a_fad.dx(i)/a_fad.val());
  COMPARE_FADS(c_fad, t2);

  c_fad = std::pow(val, b_fad);
  FadType t3(n, std::pow(val, b_fad.val()));
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) =
      std::pow(val, b_fad.val())*std::log(val)*b_fad.dx(i);
  COMPARE_FADS(c_fad, t3);

  val = 0.0;
  c_fad = std::pow(a_fad, val);
  FadType t4(n, std::pow(a_fad.val(), val));
  for (int i=0; i<n; i++)
    t4.fastAccessDx(i) = 0.0;
  COMPARE_FADS(c_fad, t4);

  c_fad = std::pow(val, b_fad);
  FadType t5(n, std::pow(val, b_fad.val()));
  for (int i=0; i<n; i++)
    t5.fastAccessDx(i) = 0.0;
  COMPARE_FADS(c_fad, t5);

  FadType aa_fad = a_fad;
  aa_fad.val() = 0.0;
  c_fad = std::pow(aa_fad, b_fad);
  FadType t6(n, std::pow(aa_fad.val(),b_fad.val()));
  for (int i=0; i<n; i++)
    t6.fastAccessDx(i) = 0.0;
  COMPARE_FADS(c_fad, t6);

  FadType bb_fad = b_fad;
  bb_fad.val() = 0.0;
  c_fad = std::pow(a_fad, bb_fad);
  FadType t7(n, std::pow(a_fad.val(),bb_fad.val()));
  for (int i=0; i<n; i++)
    t7.fastAccessDx(i) =
      std::pow(a_fad.val(),bb_fad.val())*(bb_fad.val()*a_fad.dx(i)/a_fad.val()
                                          + std::log(a_fad.val())*b_fad.dx(i));
  COMPARE_FADS(c_fad, t7);
}

TYPED_TEST_P(FadOpsUnitTest2, testEqualsLR) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_fad = a_fad;
  aa_fad = 1.0;
  aa_fad = aa_fad + b_fad;
  c_fad = 1.0 + b_fad;
  COMPARE_FADS(aa_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testPlusEqualsLR) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_fad = a_fad;
  aa_fad = 1.0;
  aa_fad += aa_fad + b_fad;
  c_fad = 1.0 + 1.0 + b_fad;
  COMPARE_FADS(aa_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testMinusEqualsLR) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_fad = a_fad;
  aa_fad = 1.0;
  aa_fad -= aa_fad + b_fad;
  c_fad = 1.0 - 1.0 - b_fad;
  COMPARE_FADS(aa_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testTimesEqualsLR) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_fad = a_fad;
  aa_fad = 2.0;
  aa_fad *= aa_fad + b_fad;
  c_fad = 2.0 * (2.0 + b_fad);
  COMPARE_FADS(aa_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testDivideEqualsLR) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_fad = a_fad;
  aa_fad = 2.0;
  aa_fad /= aa_fad + b_fad;
  c_fad = 2.0 / (2.0 + b_fad);
  COMPARE_FADS(aa_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testResizeBug6135) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;

  FadType d_fad = ScalarType(1.0);
  d_fad = d_fad + a_fad;
  c_fad = 1.0 + a_fad;
  COMPARE_FADS(d_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testEquality) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType aa_fad = a_fad;
  FadType bb_fad = b_fad;
  aa_fad.val() = 9.0;
  bb_fad.val() = 3.0;
  FadType d_fad;
  if (aa_fad == bb_fad*bb_fad)
    d_fad = aa_fad;
  c_fad = aa_fad;
  COMPARE_FADS(d_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testEqualityConstL) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType bb_fad = b_fad;
  bb_fad.val() = 3.0;
  FadType d_fad;
  if (ScalarType(9.0) == bb_fad*bb_fad)
    d_fad = a_fad;
  c_fad = a_fad;
  COMPARE_FADS(d_fad, c_fad);
}

TYPED_TEST_P(FadOpsUnitTest2, testEqualityConstR) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;

  FadType bb_fad = b_fad;
  bb_fad.val() = 3.0;
  FadType d_fad;
  if (bb_fad*bb_fad == ScalarType(9.0))
    d_fad = a_fad;
  c_fad = a_fad;
  COMPARE_FADS(d_fad, c_fad);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testLessThanOrEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  bool r1 = a_fad <= b_fad;
  bool r2 = a_fad.val() <= b_fad.val();
  ASSERT_TRUE(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = a_fad <= val;
  r2 = a_fad.val() <= val;
  ASSERT_TRUE(r1 == r2);

  r1 = val <= b_fad;
  r2 = val <= b_fad.val();
  ASSERT_TRUE(r1 == r2);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testGreaterThanOrEquals) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  bool r1 = a_fad >= b_fad;
  bool r2 = a_fad.val() >= b_fad.val();
  ASSERT_TRUE(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = a_fad >= val;
  r2 = a_fad.val() >= val;
  ASSERT_TRUE(r1 == r2);

  r1 = val >= b_fad;
  r2 = val >= b_fad.val();
  ASSERT_TRUE(r1 == r2);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testLessThan) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  bool r1 = a_fad < b_fad;
  bool r2 = a_fad.val() < b_fad.val();
  ASSERT_TRUE(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = a_fad < val;
  r2 = a_fad.val() < val;
  ASSERT_TRUE(r1 == r2);

  r1 = val < b_fad;
  r2 = val < b_fad.val();
  ASSERT_TRUE(r1 == r2);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testGreaterThan) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;

  bool r1 = a_fad > b_fad;
  bool r2 = a_fad.val() > b_fad.val();
  ASSERT_TRUE(r1 == r2);

  ScalarType val = this->urand.number();
  r1 = a_fad > val;
  r2 = a_fad.val() > val;
  ASSERT_TRUE(r1 == r2);

  r1 = val > b_fad;
  r2 = val > b_fad.val();
  ASSERT_TRUE(r1 == r2);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testACos) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::acos(a_fad);
  FadType t1(n, std::acos(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = -a_fad.dx(i)/std::sqrt(1.0 - a_fad.val()*a_fad.val());
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testASin) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::asin(a_fad);
  FadType t1(n, std::asin(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/std::sqrt(1.0 - a_fad.val()*a_fad.val());
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testATan) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::atan(a_fad);
  FadType t1(n, std::atan(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/(1.0 + a_fad.val()*a_fad.val());
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testACosh) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  FadType aa_fad = a_fad;
  if (a_fad.val() < 1.0)
    aa_fad.val() = 1.0 / a_fad.val();
  c_fad = std::acosh(aa_fad);
  FadType t1(n, std::acosh(aa_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = aa_fad.dx(i)/std::sqrt(aa_fad.val()*aa_fad.val()-1.0);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testASinh) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::asinh(a_fad);
  FadType t1(n, std::asinh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/std::sqrt(a_fad.val()*a_fad.val()+1.0);
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testATanh) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::atanh(a_fad);
  FadType t1(n, std::atanh(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = a_fad.dx(i)/(1.0 - a_fad.val()*a_fad.val());
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testAbs) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::abs(a_fad);
  FadType t1(n, std::abs(a_fad.val()));
  for (int i=0; i<n; i++) {
    if (a_fad.val() >= 0)
      t1.fastAccessDx(i) = a_fad.dx(i);
    else
      t1.fastAccessDx(i) = -a_fad.dx(i);
  }
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testFAbs) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::fabs(a_fad);
  FadType t1(n, std::fabs(a_fad.val()));
  for (int i=0; i<n; i++) {
    if (a_fad.val() >= 0)
      t1.fastAccessDx(i) = a_fad.dx(i);
    else
      t1.fastAccessDx(i) = -a_fad.dx(i);
  }
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testCbrt) {
  typedef decltype(this->a_fad_) FadType;
  auto a_fad = this->a_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::cbrt(a_fad);
  FadType t1(n, std::cbrt(a_fad.val()));
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) =
      a_fad.dx(i)/(3.*std::cbrt(a_fad.val()*a_fad.val()));
  COMPARE_FADS(c_fad, t1);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testATan2) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  c_fad = std::atan2(a_fad, b_fad);
  FadType t1(n, std::atan2(a_fad.val(),b_fad.val()));
  ScalarType t = a_fad.val()*a_fad.val() +
    b_fad.val()*b_fad.val();
  for (int i=0; i<n; i++)
    t1.fastAccessDx(i) = (b_fad.val()*a_fad.dx(i) -
                          a_fad.val()*b_fad.dx(i))/t;
  COMPARE_FADS(c_fad, t1);

  ScalarType val = this->urand.number();
  c_fad = std::atan2(a_fad, val);
  FadType t2(n, std::atan2(a_fad.val(), val));
  t = a_fad.val()*a_fad.val() + val*val;
  for (int i=0; i<n; i++)
    t2.fastAccessDx(i) = val*a_fad.dx(i)/t;
  COMPARE_FADS(c_fad, t2);

  c_fad = std::atan2(val, b_fad);
  FadType t3(n, std::atan2(val, b_fad.val()));
  t = val*val + b_fad.val()*b_fad.val();
  for (int i=0; i<n; i++)
    t3.fastAccessDx(i) = -val*b_fad.dx(i)/t;
  COMPARE_FADS(c_fad, t3);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testMax) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  ScalarType val;

  // Fad, Fad
  FadType aa_fad = a_fad + 1.0;
  c_fad = max(aa_fad, a_fad);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = max(a_fad, aa_fad);
  COMPARE_FADS(c_fad, aa_fad);

  // Expr, Fad
  c_fad = max(a_fad+1.0, a_fad);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = max(a_fad, a_fad+1.0);
  COMPARE_FADS(c_fad, aa_fad);

  // Expr, Expr (same)
  c_fad = max(a_fad+1.0, a_fad+1.0);
  COMPARE_FADS(c_fad, aa_fad);

  // Expr, Expr (different)
  c_fad = max(a_fad+1.0, a_fad-1.0);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = max(a_fad-1.0, a_fad+1.0);
  COMPARE_FADS(c_fad, aa_fad);

  // Fad, const
  val = a_fad.val() + 1;
  c_fad = max(a_fad, val);
  COMPARE_VALUES(c_fad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_fad.dx(i), 0.0);
  val = a_fad.val() - 1;
  c_fad = max(a_fad, val);
  COMPARE_FADS(c_fad, a_fad);
  val = b_fad.val() + 1;
  c_fad = max(val, b_fad);
  COMPARE_VALUES(c_fad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_fad.dx(i), 0.0);
  val = b_fad.val() - 1;
  c_fad = max(val, b_fad);
  COMPARE_FADS(c_fad, b_fad);

  // Expr, const
  val = a_fad.val();
  c_fad = max(a_fad+1.0, val);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = max(val, a_fad+1.0);
  COMPARE_FADS(c_fad, aa_fad);
}

TYPED_TEST_P(RealFadOpsUnitTest2, testMin) {
  typedef decltype(this->a_fad_) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  auto a_fad = this->a_fad_;
  auto b_fad = this->b_fad_;
  auto c_fad = this->c_fad_;
  auto n = this->n_;

  ScalarType val;

  // Fad, Fad
  FadType aa_fad = a_fad - 1.0;
  c_fad = min(aa_fad, a_fad);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = min(a_fad, aa_fad);
  COMPARE_FADS(c_fad, aa_fad);

  // Expr, Fad
  c_fad = min(a_fad-1.0, a_fad);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = min(a_fad, a_fad-1.0);
  COMPARE_FADS(c_fad, aa_fad);

  // Expr, Expr (same)
  c_fad = min(a_fad-1.0, a_fad-1.0);
  COMPARE_FADS(c_fad, aa_fad);

  // Expr, Expr (different)
  c_fad = min(a_fad+1.0, a_fad-1.0);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = min(a_fad-1.0, a_fad+1.0);
  COMPARE_FADS(c_fad, aa_fad);

  // Fad, const
  val = a_fad.val() - 1;
  c_fad = min(a_fad, val);
  COMPARE_VALUES(c_fad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_fad.dx(i), 0.0);
  val = a_fad.val() + 1;
  c_fad = min(a_fad, val);
  COMPARE_FADS(c_fad, a_fad);
  val = b_fad.val() - 1;
  c_fad = min(val, b_fad);
  COMPARE_VALUES(c_fad.val(), val);
  for (int i=0; i<n; i++)
    COMPARE_VALUES(c_fad.dx(i), 0.0);
  val = b_fad.val() + 1;
  c_fad = min(val, b_fad);
  COMPARE_FADS(c_fad, b_fad);

  // Expr, const
  val = a_fad.val();
  c_fad = min(a_fad-1.0, val);
  COMPARE_FADS(c_fad, aa_fad);
  c_fad = min(val, a_fad-1.0);
  COMPARE_FADS(c_fad, aa_fad);
}

REGISTER_TYPED_TEST_SUITE_P(
  FadOpsUnitTest2,
  testAddition,
  testSubtraction,
  testMultiplication,
  testDivision,
  testEquals,
  testNotEquals,
  testUnaryPlus,
  testUnaryMinus,
  testExp,
  testLog,
  testLog10,
  testSqrt,
  testCos,
  testSin,
  testTan,
  testCosh,
  testSinh,
  testTanh,
  testPlusEquals,
  testMinusEquals,
  testTimesEquals,
  testDivideEquals,
  testPow,
  testEqualsLR,
  testPlusEqualsLR,
  testMinusEqualsLR,
  testTimesEqualsLR,
  testDivideEqualsLR,
  testResizeBug6135,
  testEquality,
  testEqualityConstL,
  testEqualityConstR);

REGISTER_TYPED_TEST_SUITE_P(
  RealFadOpsUnitTest2,
  testLessThanOrEquals,
  testGreaterThanOrEquals,
  testLessThan,
  testGreaterThan,
  testACos,
  testASin,
  testATan,
  testACosh,
  testASinh,
  testATanh,
  testAbs,
  testFAbs,
  testCbrt,
  testATan2,
  testMax,
  testMin
  );

#endif // FADUNITTESTS2_HPP
