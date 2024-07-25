// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Random.hpp"

typedef Sacado::Fad::DFad<double> DFadType;
typedef Sacado::LFad::LogicalSparse<double,bool> LSType;

// gtest includes
#include <gtest/gtest.h>

// A class for testing each DFad operation
class LogicalSparseOpsUnitTest : public ::testing::Test {
protected:

  // DFad variables
  DFadType a_dfad, b_dfad, c_dfad;

  // Logical sparse variables
  LSType a_ls, b_ls, c_ls;

  // Random number generator
  Sacado::Random<double> urand;

  // Number of derivative components
  int n;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

  LogicalSparseOpsUnitTest() :
    urand(0.0, 1.0), n(5), tol_a(1.0e-15), tol_r(1.0e-14) {}

  void SetUp() override {
    double val;

    val = urand.number();
    a_dfad = DFadType(n,val);
    a_ls = LSType(n,val);

    val = urand.number();
    b_dfad = DFadType(n,val);
    b_ls = LSType(n,val);

    val = urand.number();
    c_dfad = val;
    c_ls = val;

    for (int i=0; i<n; i++) {
      val = urand.number();
      a_dfad.fastAccessDx(i) = val;
      a_ls.fastAccessDx(i) = 1;

      val = urand.number();
      b_dfad.fastAccessDx(i) = val;
      b_ls.fastAccessDx(i) = 1;
    }
  }

  void TearDown() override {}

  // Assert to Fad objects are the same
  void compareFads(const DFadType& x_dfad, const LSType& x_ls) {
    // Compare sizes
    ASSERT_TRUE(x_dfad.size() == x_ls.size());

    // Compare hasFastAccess
    ASSERT_TRUE(x_dfad.hasFastAccess() == x_ls.hasFastAccess());

    // Compare values
    compareDoubles(x_dfad.val(), x_ls.val());

    for (int i=0; i<x_ls.size(); i++) {

      // Compare dx
      compareDx(x_dfad.dx(i), x_ls.dx(i));

      // Compare fastAccessDx
      compareDx(x_dfad.fastAccessDx(i), x_ls.fastAccessDx(i));
    }
  }

  // Assert two doubles are the same to relative precision
  void compareDoubles(double a, double b) {
    ASSERT_TRUE( fabs(a-b) < tol_a + tol_r*fabs(a) );
  }

  // Assert two bools are the same
  void compareBools(bool a, bool b) {
    ASSERT_TRUE( a == b );
  }

  // Assert a double and bool are same (logically)
  void compareDx(double a, bool b) {
    ASSERT_TRUE( (a && b) || !(a || b) );
  }

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

}; // class LogicalSparseOpsUnitTest

#define BINARY_OP_TEST(TESTNAME,OP)            \
  TEST_F(LogicalSparseOpsUnitTest, TESTNAME) { \
    c_dfad = a_dfad OP b_dfad;                 \
    c_ls = a_ls OP b_ls;                       \
    compareFads(c_dfad, c_ls);                 \
                                               \
    double val = urand.number();               \
    c_dfad = a_dfad OP val;                    \
    c_ls = a_ls OP val;                        \
    compareFads(c_dfad, c_ls);                 \
                                               \
    c_dfad = val OP b_dfad;                    \
    c_ls = val OP b_ls;                        \
    compareFads(c_dfad, c_ls);                 \
  }

#define RELOP_TEST(TESTNAME,OP)                             \
  TEST_F(LogicalSparseOpsUnitTest, TESTNAME) {              \
    bool r1 = a_dfad OP b_dfad;                             \
    bool r2 = a_ls OP b_ls;                                 \
    ASSERT_TRUE(r1 == r2);                                  \
                                                            \
    double val = urand.number();                            \
    r1 = a_dfad OP val;                                     \
    r2 = a_ls OP val;                                       \
    ASSERT_TRUE(r1 == r2);                                  \
                                                            \
    r1 = val OP b_dfad;                                     \
    r2 = val OP b_ls;                                       \
    ASSERT_TRUE(r1 == r2);                                  \
  }

#define BINARY_FUNC_TEST(TESTNAME,FUNC)                         \
  TEST_F(LogicalSparseOpsUnitTest, TESTNAME) {                  \
    c_dfad = FUNC (a_dfad,b_dfad);                              \
    c_ls = FUNC (a_ls,b_ls);                                    \
    compareFads(c_dfad, c_ls);                                  \
                                                                \
    double val = urand.number();                                \
    c_dfad = FUNC (a_dfad,val);                                 \
    c_ls = FUNC (a_ls,val);                                     \
    compareFads(c_dfad, c_ls);                                  \
                                                                \
    c_dfad = FUNC (val,b_dfad);                                 \
    c_ls = FUNC (val,b_ls);                                     \
    compareFads(c_dfad, c_ls);                                  \
  }

#define UNARY_OP_TEST(TESTNAME,OP)                                  \
  TEST_F(LogicalSparseOpsUnitTest, TESTNAME) {                      \
    c_dfad = OP a_dfad;                                             \
    c_ls = OP a_ls;                                                 \
    compareFads(c_dfad, c_ls);                                      \
  }

#define UNARY_FUNC_TEST(TESTNAME,FUNC)                              \
  TEST_F(LogicalSparseOpsUnitTest, TESTNAME) {                      \
    c_dfad = FUNC (a_dfad);                                         \
    c_ls = FUNC (a_ls);                                             \
    compareFads(c_dfad, c_ls);                                      \
  }

#define UNARY_ASSIGNOP_TEST(TESTNAME,OP)                            \
  TEST_F(LogicalSparseOpsUnitTest, TESTNAME) {                      \
    c_dfad OP a_dfad;                                               \
    c_ls OP a_ls;                                                   \
    compareFads(c_dfad, c_ls);                                      \
                                                                    \
    double val = urand.number();                                    \
    c_dfad OP val;                                                  \
    c_ls OP val;                                                    \
    compareFads(c_dfad, c_ls);                                      \
  }

BINARY_OP_TEST(testAddition, +)
BINARY_OP_TEST(testSubtraction, -)
BINARY_OP_TEST(testMultiplication, *)
BINARY_OP_TEST(testDivision, /)

RELOP_TEST(testEquals, ==)
RELOP_TEST(testNotEquals, !=)
RELOP_TEST(testLessThanOrEquals, <=)
RELOP_TEST(testGreaterThanOrEquals, >=)
RELOP_TEST(testLessThan, <)
RELOP_TEST(testGreaterThan, >)

BINARY_FUNC_TEST(testPow, pow)

UNARY_OP_TEST(testUnaryPlus, +)
UNARY_OP_TEST(testUnaryMinus, -)

UNARY_FUNC_TEST(testExp, exp)
UNARY_FUNC_TEST(testLog, log)
UNARY_FUNC_TEST(testLog10, log10)
UNARY_FUNC_TEST(testSqrt, sqrt)
UNARY_FUNC_TEST(testCos, cos)
UNARY_FUNC_TEST(testSin, sin)
UNARY_FUNC_TEST(testTan, tan)
UNARY_FUNC_TEST(testACos, acos)
UNARY_FUNC_TEST(testASin, asin)
UNARY_FUNC_TEST(testATan, atan)
UNARY_FUNC_TEST(testCosh, cosh)
UNARY_FUNC_TEST(testSinh, sinh)
UNARY_FUNC_TEST(testTanh, tanh)
UNARY_FUNC_TEST(testAbs, abs)
UNARY_FUNC_TEST(testFAbs, fabs)

UNARY_ASSIGNOP_TEST(testPlusEquals, +=)
UNARY_ASSIGNOP_TEST(testMinusEquals, -=)
UNARY_ASSIGNOP_TEST(testTimesEquals, *=)
UNARY_ASSIGNOP_TEST(testDivideEquals, /=)

TEST_F(LogicalSparseOpsUnitTest, testComposite1) {
  c_dfad = composite1(a_dfad, b_dfad);
  c_ls = composite1(a_ls, b_ls);
  compareFads(c_dfad, c_ls);
}

TEST_F(LogicalSparseOpsUnitTest, testPlusLR) {
  DFadType aa_dfad = a_dfad;
  LSType aa_ls = a_ls;
  aa_dfad = 1.0;
  aa_ls = 1.0;
  aa_dfad = aa_dfad + b_dfad;
  aa_ls = aa_ls + b_ls;
  compareFads(aa_dfad, aa_ls);
}

TEST_F(LogicalSparseOpsUnitTest, testMinusLR) {
  DFadType aa_dfad = a_dfad;
  LSType aa_ls = a_ls;
  aa_dfad = 1.0;
  aa_ls = 1.0;
  aa_dfad = aa_dfad - b_dfad;
  aa_ls = aa_ls - b_ls;
  compareFads(aa_dfad, aa_ls);
}

TEST_F(LogicalSparseOpsUnitTest, testTimesLR) {
  DFadType aa_dfad = a_dfad;
  LSType aa_ls = a_ls;
  aa_dfad = 2.0;
  aa_ls = 2.0;
  aa_dfad = aa_dfad * b_dfad;
  aa_ls = aa_ls * b_ls;
  compareFads(aa_dfad, aa_ls);
}

TEST_F(LogicalSparseOpsUnitTest, testDivideLR) {
  DFadType aa_dfad = a_dfad;
  LSType aa_ls = a_ls;
  aa_dfad = 2.0;
  aa_ls = 2.0;
  aa_dfad = aa_dfad / b_dfad;
  aa_ls = aa_ls / b_ls;
  compareFads(aa_dfad, aa_ls);
}

TEST_F(LogicalSparseOpsUnitTest, testMax) {
  double val;

  // LFAd, LFad
  LSType aa_ls = a_ls + 1.0;
  c_ls = max(aa_ls, a_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = max(a_ls, aa_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, LFad
  c_ls = max(a_ls+1.0, a_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = max(a_ls, a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, Expr (same)
  c_ls = max(a_ls+1.0, a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, Expr (different)
  c_ls = max(a_ls+1.0, a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = max(a_ls-1.0, a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // LFad, const
  val = a_ls.val() + 1;
  c_ls = max(a_ls, val);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  val = a_ls.val() - 1;
  c_ls = max(a_ls, val);
  compareDoubles(c_ls.val(), a_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), a_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), a_ls.fastAccessDx(i));
  }
  val = b_ls.val() + 1;
  c_ls = max(val, b_ls);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  val = b_ls.val() - 1;
  c_ls = max(val, b_ls);
  compareDoubles(c_ls.val(), b_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), b_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), b_ls.fastAccessDx(i));
  }

  // Expr, const
  val = a_ls.val();
  c_ls = max(a_ls+1.0, val);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = max(val, a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
}

TEST_F(LogicalSparseOpsUnitTest, testMin) {
  double val;

  // LFad, LFad
  LSType aa_ls = a_ls - 1.0;
  c_ls = min(aa_ls, a_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = min(a_ls, aa_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, LFad
  c_ls = min(a_ls-1.0, a_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = min(a_ls, a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, Expr (same)
  c_ls = min(a_ls-1.0, a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, Expr (different)
  c_ls = min(a_ls+1.0, a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = min(a_ls-1.0, a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // LFad, const
  val = a_ls.val() - 1;
  c_ls = min(a_ls, val);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  val = a_ls.val() + 1;
  c_ls = min(a_ls, val);
  compareDoubles(c_ls.val(), a_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), a_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), a_ls.fastAccessDx(i));
  }
  val = b_ls.val() - 1;
  c_ls = min(val, b_ls);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  val = b_ls.val() + 1;
  c_ls = min(val, b_ls);
  compareDoubles(c_ls.val(), b_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), b_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), b_ls.fastAccessDx(i));
  }

  // Expr, const
  val = a_ls.val();
  c_ls = min(a_ls-1.0, val);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  c_ls = min(val, a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
}
