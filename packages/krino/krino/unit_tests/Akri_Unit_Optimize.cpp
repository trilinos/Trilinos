#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_ROLOptimize.hpp>
#include <Akri_LogRedirecter.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <vector>

#include <Akri_UnitTestUtils.hpp>

namespace krino {

using MinimizeDistributedVectorObjective =
    std::function<void(const ObjectiveInterface &objFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)>;

void optimize_using_krino_lbfgs(const ObjectiveInterface &objFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  krino::lbfgs(objFn, x, xTol, gradTol, 1.0, maxIter, 10, true);
}

void optimize_using_ROL_lbfgs(const ObjectiveInterface &objFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  rol_optimize(objFn, x, xTol, gradTol, 1.0, maxIter, true);
}

std::tuple<double,double> minimize_function_along_gradient(const ObjectiveInterface &objFn,
    const MinimizeDistributedVectorObjective & minimize,
    const double x0,
    const int maxIter)
{
  DistributedVector soln{x0};
  std::cout << "Minimization starting at " << x0 << " " << objFn.compute_value(soln) << std::endl;
  const double xTol = 1.e-6;
  const double gradTol = 1.e-6;
  minimize(objFn, soln, xTol, gradTol, maxIter);
  double fmin = objFn.compute_value(soln);
  std::cout << "soln " << soln[0] << " " << fmin << std::endl;
  return {soln[0], fmin};
}

void expect_find_minimum_along_gradient(const ObjectiveInterface &objFn,
    const MinimizeDistributedVectorObjective & minimize,
    const double x0,
    const double goldXmin,
    const double xtol,
    const double ftol,
    const int maxIter)
{
  const auto & [xmin, fmin] = minimize_function_along_gradient(objFn, minimize, x0, maxIter);
  EXPECT_NEAR(goldXmin, xmin, xtol);
  EXPECT_LT(fmin, ftol);
}

class ScalarObjective : public ObjectiveInterface
{
public:
  virtual double compute_value(const DistributedVector & x) const override { return value(x[0]); }
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override { g = {gradient(x[0])}; }
  virtual double value(const double x) const = 0;
  virtual double gradient(const double x) const = 0;
};

class QuadraticFunction : public ScalarObjective
{
public:
  virtual double value(const double x) const override { return std::pow(x-5.,2); }
  virtual double gradient(const double x) const override { return 2*(x-5.); }
};

TEST(minimize_function, quadraticFunction)
{
  const int maxIter = 2;
  const double xtol = 1.e-6;
  const double ftol = 1.e-6;
  const QuadraticFunction objFn;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(objFn, minimizer, -5, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 0., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 10., 5., xtol, ftol, maxIter);
  }
}

class AbsFunction : public ScalarObjective
{
public:
  virtual double value(const double x) const override { return std::pow(std::abs(x-5.)/5., 3); }
  virtual double gradient(const double x) const override
  { return (x<5.) ? (-3/5.*std::pow((5.-x)/5.,2)) : (3/5.*std::pow((x-5.)/5.,2)); }
};

TEST(minimize_function, absFunction)
{
  const int maxIter = 20;
  const double xtol = 1.e-2; // larger xtol because function is flatter near minimum
  const double ftol = 1.e-6;
  const AbsFunction objFn;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(objFn, minimizer, 1.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 2.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 8.1, 5., xtol, ftol, maxIter);
  }
}

class QuarticFunction : public ScalarObjective
{
public:
  virtual double value(const double x) const override { return std::pow((x-5.)/5., 4); }
  virtual double gradient(const double x) const override { return 0.8*pow((x-5.)/5.,3); }
};

TEST(minimize_function, quarticFunction)
{
  const int maxIter = 20;
  const double xtol = 1.e-1; // large xtol because function is so flat near minimum
  const double ftol = 1.e-6;
  const QuarticFunction objFn;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(objFn, minimizer, 1.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 2.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 8.1, 5., xtol, ftol, maxIter);
  }
}

class RosenbrockFunction : public ObjectiveInterface
{
public:
  static constexpr double alpha = 100.;
  virtual double compute_value(const DistributedVector & x) const override
  {
    double val = 0;
    for( unsigned i=0; i<x.size()/2; i++ )
    {
      val += alpha * std::pow(std::pow(x[2*i],2) - x[2*i+1], 2);
      val += std::pow(x[2*i] - 1.0, 2);
    }
    return val;
  }
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override
  {
    g.resize(x.sizes());
    for( unsigned i=0; i<x.size()/2; i++ )
    {
      g[2*i]   =  4.0*alpha*(std::pow(x[2*i],2) - x[2*i+1])*x[2*i] + 2.0*(x[2*i]-1.0);
      g[2*i+1] = -2.0*alpha*(std::pow(x[2*i],2) - x[2*i+1]);
    }
  }
  DistributedVector get_initial_guess(const unsigned n) const
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
};

void test_rosenbrock(const unsigned n, const MinimizeDistributedVectorObjective & mininimize)
{
  const double xTol = 1.e-6;
  const double gradTol = 1.e-12;
  const unsigned maxIter = 50;
  const RosenbrockFunction objFn;
  auto x = objFn.get_initial_guess(n);
  mininimize(objFn, x, xTol, gradTol, maxIter);
  const std::vector<double> gold(n, 1.0);
  expect_near_absolute(gold, x.get(), xTol);
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



