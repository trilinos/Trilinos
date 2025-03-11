#include "stk_middle_mesh/bspline_approximation_1d.hpp"
#include "gtest/gtest.h"
#include <fstream>

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

std::vector<double> get_x_uniform(double xmin, double xmax, int npts)
{
  std::vector<double> xVals(npts);
  double deltaX = (xmax - xmin) / (npts - 1);

  for (int i = 0; i < npts; ++i)
    xVals[i] = xmin + deltaX * i;

  return xVals;
}

template <typename T>
std::vector<double> get_f_values(const std::vector<double>& xValues, T func)
{
  std::vector<double> fVals(xValues.size());
  for (size_t i = 0; i < xValues.size(); ++i)
    fVals[i] = func(xValues[i]);

  return fVals;
}

} // namespace

TEST(BSplineApproximation1D, Zero)
{
  auto xVals = get_x_uniform(0, 1, 4);
  auto fVals = get_f_values(xVals, [](double /*x*/) { return 0; });

  utils::impl::BSplineApproximation1D approximator(4);
  approximator.approximate(xVals, fVals);

  auto xTest = get_x_uniform(0, 1, 100);
  for (auto x : xTest)
    EXPECT_NEAR(approximator.eval(x), 0, 1e-13);
}

TEST(BSplineApproximation1D, Constant)
{
  auto xVals = get_x_uniform(0, 1, 20);
  // std::vector<double> x_vals = {0, 1.0/3.0, 2.0/3.0, 1};
  auto fVals = get_f_values(xVals, [](double /*x*/) { return 2; });

  utils::impl::BSplineApproximation1D approximator(4);
  approximator.approximate(xVals, fVals);

  auto xTest = get_x_uniform(0, 1, 100);
  std::vector<double> fValsApprox;
  for (auto x : xTest)
  {
    fValsApprox.push_back(approximator.eval(x));
    EXPECT_NEAR(approximator.eval(x), 2, 1e-13);
  }
}

TEST(BSplineApproximation1D, Linear)
{
  auto xVals = get_x_uniform(0, 1, 100);
  // std::vector<double> x_vals = {0, 1.0/3.0, 2.0/3.0, 1};
  auto fVals = get_f_values(xVals, [](double x) { return 2 * x; });

  utils::impl::BSplineApproximation1D approximator(40);
  approximator.approximate(xVals, fVals);

  auto xTest = get_x_uniform(0, 1, 100);
  std::vector<double> fValsApprox;
  for (auto x : xTest)
  {
    fValsApprox.push_back(approximator.eval(x));
    EXPECT_NEAR(approximator.eval(x), 2 * x, 1e-13);
  }
}

TEST(BSplineApproximation1D, Quadratic)
{
  auto xVals = get_x_uniform(0, 1, 100);
  // std::vector<double> x_vals = {0, 1.0/3.0, 2.0/3.0, 1};
  auto fVals = get_f_values(xVals, [](double x) { return 2 * x * x; });

  utils::impl::BSplineApproximation1D approximator(40);
  approximator.approximate(xVals, fVals);

  auto xTest = get_x_uniform(0, 1, 100);
  std::vector<double> fValsApprox;
  for (auto x : xTest)
  {
    fValsApprox.push_back(approximator.eval(x));
    EXPECT_NEAR(approximator.eval(x), 2 * x * x, 1e-13);
  }
}

TEST(BSplineApproximation1D, Cubic)
{
  auto xVals = get_x_uniform(0, 1, 20);
  // std::vector<double> x_vals = {0, 1.0/3.0, 2.0/3.0, 1};
  auto fVals = get_f_values(xVals, [](double x) { return 2 * x * x * x; });

  utils::impl::BSplineApproximation1D approximator(4);
  approximator.approximate(xVals, fVals);

  auto xTest = get_x_uniform(0, 1, 100);
  std::vector<double> fValsApprox;
  for (auto x : xTest)
  {
    fValsApprox.push_back(approximator.eval(x));
    EXPECT_NEAR(approximator.eval(x), 2 * x * x * x, 1e-13);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
