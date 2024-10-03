// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// gtest includes
#include <gtest/gtest.h>

struct CompareFloats {
public:
  double tol_a, tol_r;

  CompareFloats(double tol_a_, double tol_r_) : tol_a(tol_a_), tol_r(tol_r_) {}

  template <typename ScalarType>
  bool operator() (const ScalarType& a, const ScalarType& b) {
    return std::abs(a-b) < tol_a + tol_r*std::abs(a);
  }
};

struct CompareFads {
public:
  CompareFloats cmp;

  CompareFads(double tol_a, double tol_r) : cmp(tol_a, tol_r) {}

  template <typename FadType1, typename FadType2>
  bool operator() (const FadType1& a, const FadType2& b)
  {
    if (a.size() != b.size()) return false;
    if (a.hasFastAccess() != b.hasFastAccess()) return false;
    if (!cmp(a.val(), b.val())) return false;
    for (int i=0; i<a.size(); ++i) {
      if (!cmp(a.dx(i), b.dx(i))) return false;
      if (!cmp(a.fastAccessDx(i), b.fastAccessDx(i))) return false;
    }
    return true;
  }

  template <typename FadType1, typename FadType2>
  ::testing::AssertionResult operator() (const char* a_expr, const char* b_expr,
                                         const FadType1& a, const FadType2& b)
  {
    bool success = (*this)(a,b);
    if (success) return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure()
      << "Fad's do not match!" << std::endl
      << a_expr << " = " << a << std::endl
      << b_expr << " = " << b << std::endl;
  }
};

struct CompareNestedFads {
public:
  CompareFads cmp;

  CompareNestedFads(double tol_a, double tol_r) : cmp(tol_a, tol_r) {}

  template <typename FadType1, typename FadType2>
  bool operator() (const FadType1& a, const FadType2& b)
  {
    if (a.size() != b.size()) return false;
    if (a.hasFastAccess() != b.hasFastAccess()) return false;
    if (!cmp(a.val(), b.val())) return false;
    for (int i=0; i<a.size(); ++i) {
      if (!cmp(a.dx(i), b.dx(i))) return false;
      if (!cmp(a.fastAccessDx(i), b.fastAccessDx(i))) return false;
    }
    return true;
  }

  template <typename FadType1, typename FadType2>
  ::testing::AssertionResult operator() (const char* a_expr, const char* b_expr,
                                         const FadType1& a, const FadType2& b)
  {
    bool success = (*this)(a,b);
    if (success) return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure()
      << "Fad's do not match!" << std::endl
      << a_expr << " = " << a << std::endl
      << b_expr << " = " << b << std::endl;
  }
};

#define COMPARE_VALUES(a, b)                                    \
  ASSERT_PRED2(CompareFloats(this->tol_a, this->tol_r), a, b);

#define COMPARE_FADS(a, b)                                      \
  ASSERT_PRED_FORMAT2(CompareFads(this->tol_a, this->tol_r), a, b);

#define COMPARE_NESTED_FADS(a, b)                               \
  ASSERT_PRED_FORMAT2(CompareNestedFads(this->tol_a, this->tol_r), a, b);
