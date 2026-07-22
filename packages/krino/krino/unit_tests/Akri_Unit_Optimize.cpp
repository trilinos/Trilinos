#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_ROLOptimize.hpp>
#include <Akri_Unit_LogRedirecter.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <vector>

#include <Akri_UnitTestUtils.hpp>

namespace krino {

using MinimizeDistributedVectorObjective =
    std::function<void(const DistributedVectorObjectiveFn &objFn,
        const DistributedVectorObjectiveSensFn &gradObjFn,
        DistributedVector& x,
        const double xTol,
        const double gradTol,
        const unsigned maxIter)>;

void optimize_using_krino_lbfgs(const DistributedVectorObjectiveFn &objFn,
    const DistributedVectorObjectiveSensFn &gradObjFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  krino::lbfgs(objFn, gradObjFn, x, xTol, gradTol, maxIter);
}

void optimize_using_ROL_lbfgs(const DistributedVectorObjectiveFn &objFn,
    const DistributedVectorObjectiveSensFn &gradObjFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  rol_optimize(objFn, gradObjFn, x, xTol, gradTol, maxIter);
}

std::tuple<double,double> minimize_function_along_gradient(const std::function<double(double)> & fn,
    const std::function<double(double)> & dfdx,
    const MinimizeDistributedVectorObjective & minimize,
    const double x0,
    const int maxIter)
{
  const auto vecfn = [&](const DistributedVector& x) { return fn(x[0]); };
  const auto vecgrad = [&](const DistributedVector& x, DistributedVector& grad) { grad = {dfdx(x[0])}; };
  std::cout << "Minimization starting at " << x0 << " " << fn(x0) << std::endl;
  const double xTol = 1.e-6;
  const double gradTol = 1.e-6;
  DistributedVector soln{x0};
  minimize(vecfn, vecgrad, soln, xTol, gradTol, maxIter);
  double xmin = soln[0];
  double fmin = fn(xmin);
  std::cout << "soln " << soln[0] << " " << fmin << std::endl;
  return {xmin, fmin};
}

void expect_find_minimum_along_gradient(const std::function<double(double)> & fn,
    const std::function<double(double)> & dfdx,
    const MinimizeDistributedVectorObjective & minimize,
    const double x0,
    const double goldXmin,
    const double xtol,
    const double ftol,
    const int maxIter)
{
  const auto & [xmin, fmin] = minimize_function_along_gradient(fn, dfdx, minimize, x0, maxIter);
  EXPECT_NEAR(goldXmin, xmin, xtol);
  EXPECT_LT(fmin, ftol);
}

TEST(minimize_function, quadraticFunction)
{
  const int maxIter = 2;
  const double xtol = 1.e-6;
  const double ftol = 1.e-6;
  const auto f = [](const double x) { return std::pow(x-5.,2); };
  const auto dfdx = [](const double x) { return 2.*(x-5.); };

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(f, dfdx, minimizer, -5, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 0., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 10., 5., xtol, ftol, maxIter);
  }
}

TEST(minimize_function, absFunction)
{
  const int maxIter = 20;
  const double xtol = 1.e-2; // larger xtol because function is flatter near minimum
  const double ftol = 1.e-6;
  const auto f = [](const double x) { return std::pow(std::abs(x-5.)/5., 3); };
  const auto dfdx = [](const double x) { return (x<5.) ? (-3/5.*std::pow((5.-x)/5.,2)) : (3/5.*std::pow((x-5.)/5.,2)); };

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 1.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 2.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 8.1, 5., xtol, ftol, maxIter);
  }
}

TEST(minimize_function, quarticFunction)
{
  const int maxIter = 20;
  const double xtol = 1.e-1; // large xtol because function is so flat near minimum
  const double ftol = 1.e-6;
  const auto f = [](const double x) { return std::pow((x-5.)/5., 4); };
  const auto dfdx = [](const double x) { return 0.8*pow((x-5.)/5.,3); };

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 1.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 2.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(f, dfdx, minimizer, 8.1, 5., xtol, ftol, maxIter);
  }
}

namespace Rosenbrock
{
double alpha = 100.;

double compute_value(const DistributedVector &x)
{
  double val = 0;
  for( unsigned i=0; i<x.size()/2; i++ )
  {
    val += alpha * std::pow(std::pow(x[2*i],2) - x[2*i+1], 2);
    val += std::pow(x[2*i] - 1.0, 2);
  }
  return val;
}

void fill_gradient(const DistributedVector &x, DistributedVector &g)
{
  g.resize(x.sizes());
  for( unsigned i=0; i<x.size()/2; i++ )
  {
    g[2*i]   =  4.0*alpha*(std::pow(x[2*i],2) - x[2*i+1])*x[2*i] + 2.0*(x[2*i]-1.0);
    g[2*i+1] = -2.0*alpha*(std::pow(x[2*i],2) - x[2*i+1]);
  }
}

DistributedVector get_initial_guess(const unsigned n)
{
  // Get Initial Guess
  DistributedVector x0(n);
  for( unsigned i=0; i<n/2; i++ )
  {
    x0[2*i]   =  1.2;
    x0[2*i+1] = -1.0;
  }
  return x0;
}
}

void test_rosenbrock(const unsigned n, const MinimizeDistributedVectorObjective & mininimize)
{
  const double xTol = 1.e-6;
  const double gradTol = 1.e-6;
  const unsigned maxIter = 50;
  auto x = Rosenbrock::get_initial_guess(n);
  mininimize(Rosenbrock::compute_value, Rosenbrock::fill_gradient, x, xTol, gradTol, maxIter);
  const std::vector<double> gold(n, 1.0);
  expect_near_absolute(gold, x.get());
}

TEST(minimize_function, rosenbrockFunction)
{
  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    test_rosenbrock(2, minimizer);
    test_rosenbrock(40, minimizer);
  }
}

}



