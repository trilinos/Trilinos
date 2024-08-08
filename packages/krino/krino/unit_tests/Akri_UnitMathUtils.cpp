// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_DiagWriter.hpp>
#include <Akri_MathUtil.hpp>
#include <gtest/gtest.h>
#include <functional>

namespace krino {

void expect_root(const double goldRoot, const double xTol, const std::function<double(const double)> & f)
{
  const unsigned maxIters = 100;
  const auto result = find_root(f, 0., 1., f(0.), f(1.), maxIters, xTol);
  ASSERT_TRUE(result.first);
  EXPECT_NEAR(goldRoot, result.second, xTol);
}

void expect_root_newton_raphson(const double goldRoot, const double guess, const double fTol, const std::function<std::pair<double,double>(const double)> & f)
{
  const unsigned maxIters = 100;
  const auto result = find_root_newton_raphson(f, guess, maxIters, fTol);
  ASSERT_TRUE(result.first);
  const auto valueAndDeriv = f(result.second);
  const double xTol = fTol / std::abs(valueAndDeriv.second);
  EXPECT_NEAR(0., valueAndDeriv.first, fTol);
  EXPECT_NEAR(goldRoot, result.second, xTol);
}

TEST(find_root, givenPolynomialFunction_findRootWithinTolerance)
{
  const double tol = 1.e-5;
  expect_root(0.25, tol, [](const double x){ return x-0.25; });
  expect_root(0.25, tol, [](const double x){ return x*x-0.25*0.25; });
  expect_root(0.25, tol, [](const double x){ return x*x*x-0.25*0.25*0.25; });
}


TEST(find_root_newton_raphson, givenPolynomialFunctionWithCorrectJacobian_findRootWithinTolerance)
{
  const double tol = 1.e-5;
  const double guess = 1.;
  expect_root_newton_raphson(0.25, guess, tol, [](const double x){ std::cout << "Eval at " << x << std::endl; return std::make_pair(x-0.25, 1.); });
  expect_root_newton_raphson(0.25, guess, tol, [](const double x){ std::cout << "Eval at " << x << std::endl; return std::make_pair(x*x-0.25*0.25, 2.*x); });
  expect_root_newton_raphson(0.25, guess, tol, [](const double x){ std::cout << "Eval at " << x << std::endl; return std::make_pair(x*x*x-0.25*0.25*0.25, 3.*x*x); });
}

TEST(find_root_newton_raphson, givenPolynomialFunctionWithWRONGJacobian_findRootWithinTolerance)
{
  const double tol = 1.e-5;
  const double guess = 1.;
  const double error = 0.1; // Less than 1 to make function overshoot
  expect_root_newton_raphson(0.25, guess, tol, [error](const double x){ std::cout << "Eval at " << x << std::endl; return std::make_pair(x-0.25, 1.*error); });
  expect_root_newton_raphson(0.25, guess, tol, [error](const double x){ std::cout << "Eval at " << x << std::endl; return std::make_pair(x*x-0.25*0.25, 2.*x*error); });
  expect_root_newton_raphson(0.25, guess, tol, [error](const double x){ std::cout << "Eval at " << x << std::endl; return std::make_pair(x*x*x-0.25*0.25*0.25, 3.*x*x*error); });
}

void expect_quadratic_crossing(const double gold, const std::array<double,3> & edgeVals)
{
  EXPECT_NEAR(gold, find_quadratic_crossing(edgeVals[0],edgeVals[1],edgeVals[2]), 1.e-6);
}

std::array<double,3> compute_edge_values(const double crossing1, const double crossing2)
{
  std::array<double,3> edgeVals = {{(0.-crossing1)*(0.-crossing2), (1.-crossing1)*(1.-crossing2), (0.5-crossing1)*(0.5-crossing2)}};
  return edgeVals;
}

TEST(compute_edge_values, singleCrossings)
{
  expect_quadratic_crossing(0.25, compute_edge_values(0.25, -0.25));
  expect_quadratic_crossing(0.25, compute_edge_values(0.25, 1.25));
  expect_quadratic_crossing(0.33, compute_edge_values(0.33, -0.25));
  expect_quadratic_crossing(0.77, compute_edge_values(0.77, 1.25));
}

}

