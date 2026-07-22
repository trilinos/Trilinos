// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_UNIT_TESTS_INCLUDE_AKRI_UNITTESTUTILS_H_
#define KRINO_UNIT_TESTS_INCLUDE_AKRI_UNITTESTUTILS_H_
#include <vector>
#include <stk_math/StkVector.hpp>
#include <gtest/gtest.h>

namespace krino {

class Facet2d;
class Facet3d;

void expect_eq(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double relativeTol=1.e-6);
void expect_eq_absolute(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double absoluteTol=1.e-6);
inline void expect_near(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double relativeTol=1.e-6) { return expect_eq(gold, result, relativeTol); }
inline void expect_near_absolute(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double absoluteTol=1.e-6) { return expect_eq_absolute(gold, result, absoluteTol); }
void expect_eq(const Facet2d & gold, const Facet2d & result, const double relativeTol=1.e-6);
void expect_eq(const Facet3d & gold, const Facet3d & result, const double relativeTol=1.e-6);
bool is_near_relative(const Facet2d & gold, const Facet2d & result, const double relativeTol=1.e-6);
bool is_near_relative(const Facet3d & gold, const Facet3d & result, const double relativeTol=1.e-6);
void expect_near_relative(const std::vector<double> & gold, const std::vector<double> & result, const double relativeTol=1.e-6);
void expect_near_absolute(const std::vector<double> & gold, const std::vector<double> & result, const double absoluteTol=1.e-6);

bool is_debug();

int num_random_test_cases(const int numDebugCases, const int numOptimizedCases);

class OptimizedOnlyTest : public ::testing::Test
{
protected:
void SetUp() override
{
  if (is_debug())
  {
    GTEST_SKIP() << "Skipping all debug build tests for this fixture";
  }
}
};

#define SkipTestMsgIf(expr, message)    \
  do {                                  \
    if (expr)                           \
    {                                   \
      GTEST_SKIP() << message;          \
    }                                   \
  } while (false)

#define SkipTestIf(expr) SkipTestMsgIf(expr, "")

}

#endif /* KRINO_UNIT_TESTS_INCLUDE_AKRI_UNITTESTUTILS_H_ */
