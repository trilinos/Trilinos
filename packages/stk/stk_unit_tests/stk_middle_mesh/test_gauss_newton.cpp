#include "gtest/gtest.h"
#include "stk_middle_mesh/gauss_newton.hpp"

using namespace stk::middle_mesh::utils::impl;

TEST(GaussNewton, Quadratic)
{

  auto f = [](const std::vector<double>& x, std::vector<double>& residuals)
  {
    for (int i=0; i < 2; ++i)
      residuals[i] = x[i]*x[i];
  };

  auto jacTest = [](const std::vector<double>& x, Matrix<double>& jac)
  {
    jac.fill(0);
    for (int i=0; i < 2; ++i)
      jac(i, i) = 2*x[i];
  };

  GaussNewton gaussNewton(2, 2, 1e-40, 100);

  std::vector<double> x0 = {1, 2};
  gaussNewton.solve(f, jacTest, x0);

  EXPECT_NEAR(x0[0], 0, 1e-12);
  EXPECT_NEAR(x0[1], 0, 1e-12);
}

TEST(GaussNewton, RectangularDegenerate)
{

  auto f = [](const std::vector<double>& x, std::vector<double>& residuals)
  {
    residuals[0] = x[0] - 0.5;
    residuals[1] = x[1] - 1;
    residuals[2] = 1.5 - x[0] - x[1];
  };

  auto jacTest = [](const std::vector<double>& /*x*/, Matrix<double>& jac)
  {
    jac.fill(0);
    jac(0, 0) = 1;
    jac(1, 1) = 1;
    jac(2, 0) = -1;  jac(2, 1) = -1;
  };

  GaussNewton gaussNewton(3, 2, 1e-40, 100);

  std::vector<double> x0 = {1, 2};
  gaussNewton.solve(f, jacTest, x0);

  EXPECT_NEAR(x0[0], 0.5, 1e-12);
  EXPECT_NEAR(x0[1], 1.0, 1e-12);
}