#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_ROLOptimize.hpp>
#include <Akri_LogRedirecter.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/ParallelReduce.hpp>
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

using MinimizeVector3dObjective =
    std::function<void(const Objective3DInterface &objFn,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)>;

void optimize_vector3d_using_krino_lbfgs(const Objective3DInterface &objFn,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  krino::lbfgs(objFn, x, xTol, gradTol, 1.0, maxIter, 10, true);
}

void optimize_vector3d_using_ROL_lbfgs(const Objective3DInterface &objFn,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  rol_optimize(objFn, x, xTol, gradTol, 1.0, maxIter, true);
}

// Anisotropic quadratic bowl f(x) = sum_i w_i*(x_i - c_i)^2, minimized at c.
class Vector3dQuadraticObjective : public Objective3DInterface
{
public:
  Vector3dQuadraticObjective(const stk::math::Vector3d & center, const stk::math::Vector3d & weights)
  : myCenter(center), myWeights(weights) {}

  virtual double compute_value(const stk::math::Vector3d & x) const override
  {
    double val = 0.;
    for (unsigned i=0; i<3; ++i)
      val += myWeights[i] * std::pow(x[i] - myCenter[i], 2);
    return val;
  }
  virtual void fill_gradient(const stk::math::Vector3d & x, stk::math::Vector3d & g) const override
  {
    for (unsigned i=0; i<3; ++i)
      g[i] = 2. * myWeights[i] * (x[i] - myCenter[i]);
  }
private:
  stk::math::Vector3d myCenter;
  stk::math::Vector3d myWeights;
};

TEST(minimize_function, vector3dQuadraticFunction)
{
  const double xTol = 1.e-8;
  const double gradTol = 1.e-10;
  const unsigned maxIter = 50;
  const stk::math::Vector3d center(1.5, -2.0, 3.25);
  const Vector3dQuadraticObjective objFn(center, stk::math::Vector3d(1.0, 4.0, 0.5));

  for (auto minimizer : {optimize_vector3d_using_krino_lbfgs, optimize_vector3d_using_ROL_lbfgs})
  {
    stk::math::Vector3d x(-3.0, 5.0, 0.0);
    minimizer(objFn, x, xTol, gradTol, maxIter);
    EXPECT_NEAR(center[0], x[0], 1.e-5);
    EXPECT_NEAR(center[1], x[1], 1.e-5);
    EXPECT_NEAR(center[2], x[2], 1.e-5);
    EXPECT_LT(objFn.compute_value(x), 1.e-8);
  }
}

// Separable quadratic with one DOF shared (ghosted) across two ranks:
//   f(x) = (x0-1)^2 + (x1-2)^2 + 2*(xS-3)^2
// Rank 0 owns x0 and the shared DOF xS; rank 1 owns x1 and holds a ghost of xS.
// Each rank contributes one (xS-3)^2 term, so the assembled gradient of xS sums
// contributions from both ranks (mirrored back onto the ghost) -- exactly the
// owned/ghost invariant that Dot/xpby preserve. This distinguishes the
// parallel-consistent ROLDistributedVector path (correct) from a plain
// StdVector, whose un-reduced dot/norm would diverge between ranks.
class SeparableSharedObjective : public ObjectiveInterface
{
public:
  virtual double compute_value(const DistributedVector & x) const override
  {
    const int rank = stk::parallel_machine_rank(x.comm());
    const double target = (rank == 0) ? 1.0 : 2.0; // target for the owned non-shared DOF (index 0)
    const double localVal = std::pow(x[0] - target, 2) + std::pow(x[1] - 3.0, 2); // index 1 is xS (owned on r0, ghost on r1)
    double globalVal = 0.0;
    stk::all_reduce_sum(x.comm(), &localVal, &globalVal, 1);
    return globalVal;
  }

  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override
  {
    const int rank = stk::parallel_machine_rank(x.comm());
    g.assign(x.comm(), x.size(), x.local_size(), 0.0);
    const double target = (rank == 0) ? 1.0 : 2.0;
    g[0] = 2.0 * (x[0] - target);
    // Assemble the shared-DOF gradient: sum each rank's contribution, then mirror
    // the result onto both the owner (r0 index 1) and the ghost (r1 index 1).
    const double localSharedContrib = 2.0 * (x[1] - 3.0);
    double globalSharedGrad = 0.0;
    stk::all_reduce_sum(x.comm(), &localSharedContrib, &globalSharedGrad, 1);
    g[1] = globalSharedGrad;
  }
};

TEST(minimize_function, separableSharedDistributedVectorNp2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) GTEST_SKIP();
  const int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

  const SeparableSharedObjective objFn;

  const double xTol = 1.e-8;
  const double gradTol = 1.e-10;
  const unsigned maxIter = 50;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    // rank 0: owned {x0, xS}; rank 1: owned {x1}, ghost {xS}. Ghost starts consistent.
    DistributedVector x = (rank == 0)
        ? DistributedVector(MPI_COMM_WORLD, std::vector<double>{0.0, 0.0}, std::vector<double>{})
        : DistributedVector(MPI_COMM_WORLD, std::vector<double>{0.0}, std::vector<double>{0.0});

    minimizer(objFn, x, xTol, gradTol, maxIter);

    if (rank == 0)
    {
      EXPECT_NEAR(1.0, x[0], 1.e-5); // owned x0
      EXPECT_NEAR(3.0, x[1], 1.e-5); // owned shared xS
    }
    else
    {
      EXPECT_NEAR(2.0, x[0], 1.e-5); // owned x1
      EXPECT_NEAR(3.0, x[1], 1.e-5); // ghost xS mirrors owner
    }
    EXPECT_LT(objFn.compute_value(x), 1.e-8);
  }
}

}



